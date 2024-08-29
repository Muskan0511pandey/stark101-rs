pub mod channel;
pub mod field;
pub mod merkle_tree;
pub mod polynomial;
pub mod utils;
pub mod verifier;

use channel::Channel;
use chrono::Local;
use field::{Field, FieldElement};
use log::{Level, LevelFilter, Metadata, Record};
use merkle_tree::MerkleTree;
use polynomial::{interpolate_lagrange_polynomials, Polynomial};
use utils::{decommit_fri, fri_commit, generate_eval_domain_for_trace, generate_trace};
use verifier::verify_proof;
static CONSOLE_LOGGER: ConsoleLogger = ConsoleLogger;
struct ConsoleLogger;

impl log::Log for ConsoleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Debug
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!(
                "{} [{}] {}:{} - {}",
                Local::now().format("%Y-%m-%dT%H:%M:%S"),
                record.level(),
                record.module_path().unwrap(),
                record.line().unwrap(),
                record.args()
            );
        }
    }

    fn flush(&self) {}
}
// Stark 101 from stark ware written in rust.
// Fibonacci Sq Mod Prime
// a{n+2} = a{n+1}^2 + a{n}^2 mod prime
// prime = 3.2^30 + 1 = 3221225473

// Fibonacci Sq Mod 3221225473 with a{0} = 1, a{1} = x, we have a{1022} = 2338775057
fn main() {
    log::set_logger(&CONSOLE_LOGGER).unwrap();
    log::set_max_level(LevelFilter::Info);

    // Part 1: LDE and Commitment
    // LDE: Low Degree Extension in 3 steps.
    // 1. Generate Input.
    // 2. Interpolate.
    // 3. Extend i.e Evalute at many points.
    // LDE for STARK: Input a0, a1, a2, ..., a1022
    //                Evaluation domain: 1, g, g^2 g^3 ... g^1022
    //                g - element from F_p
    let prime_mod =3221225473;
    let field =Field::new(prime_mod);
    let start_time =Local::now();
    log::info!("generating trace");
    let a0=FieldElement::new(1,field);
    let a1= FieldElement::new(3141592,field);
    let trace =generate_trace(a0,a1,1023);
    let g=FieldElement::new(5,field).pow(3*2u64.pow(20));
    let eval_domain=generate_eval_domain_for_trace(&trace, g);
    // interpoltae
    log::info!("interpolating");
    let trace_polynomial =interpolate_lagrange_polynomials(eval_domain[..1023].to_vec(),trace);
    //Extending
    log::info!("extending");
    let w=FieldElement::new(5,field);
    let two =FieldElement(2,field);
    let exp=two.pow(30)*FieldElement(3,field);
    let h= w.pow(exp.0);
    #[allow(non_snake_case)]
    let H: Vec<FieldElement> = (0..8192).into_iter().map(|i| h.pow(i)).collect();
    let eval_domain:Vec<FieldElement>=H.into_iter().map(|h|w*h).collect();
    let f_evaluations:Vec<FieldElement>=eval_domain
    .clone()
    .into_iter()
    .map(|h|trace_polynomial.evaluate(h
        
    ))
    .collect();
//commiting to LDE
log::info!("commiting to LDE");
let f_merkle_tree=MerkleTree::new(&f_evaluations);
let mut channel = Channel::new();
//sending merkle root to channel i.e. verifier for converting interactive to non-interactive  proof we have to use fiat shemiar heuristic
channel.send(f_merkle_tree.inner.root().unwrap().to_vec());
//generating polynomial commitments
log::info!("generating constraint polynomials");
log::debug!("generating p_0");
// p_0 = (f(x) - 1) / (x - g^0)
let p_0=(trace_polynomial.clone())-Polynomial::new_from_coefficients(vec![FieldElement(1,field)])/Polynomial::new_from_coefficients(vec![-g.pow(0),FieldElement(1,field)]);
// p_1 = (f(x) - 2338775057) / (x - g^1022)
let p_1=(trace_polynomial.clone())-Polynomial::new_from_coefficients(vec![FieldElement(2338775057,field)])/Polynomial::new_from_coefficients(vec![-g.pow(1022),FieldElement(1,field)]);
let mut  x_1024=Polynomial::new_from_coefficients(vec![FieldElement::new(0,field);1025]);
x_1024.coefficients[0]=-FieldElement::new(1,field);
x_1024.coefficients[1024]=FieldElement::new(1,field);
// p_2 = (f(x) - cp(x)) / (x - g^1023)
let p_2_den = x_1024/Polynomial::new_from_coefficients(vec![-g.pow(1021),FieldElement(1,field)])* Polynomial::new_from_coefficients(vec![-g.pow(1022),FieldElement(1,field)])*Polynomial::new_from_coefficients(vec![-g.pow(1023),FieldElement(1,field)]);
let trace_poly_g2=trace_polynomial.clone().compose(g.pow(2));
let trace_poly_g= trace_polynomial.clone().compose(g);
let square_trace_poly=trace_polynomial.clone()*trace_polynomial.clone();
let p2_num=trace_poly_g2-trace_poly_g-square_trace_poly;
let p_2= p2_num/p_2_den ;
// taking random field elements alpha_0,alpha_1,alpha_2  from channel i.e. verifier
log::info!("receive alpha_0,aplha_1,alpha_2");
let alpha_0=channel.receive_random_field_element(field);
let alpha_1=channel.receive_random_field_element(field);
let alpha_2=channel.receive_random_field_element(field);
//constructing composite polynomial cp = alpha_0*p_0 + alpha_1*p_1 + alpha_2*p_2
let cp = Polynomial::new_from_coefficients(vec![alpha_0])*p_0 + Polynomial::new_from_coefficients(vec![alpha_1])*p_1 + Polynomial::new_from_coefficients(vec![alpha_2])*p_2;
//commiting to constraint polynomial
let cp_evaluations:Vec<FieldElement>=eval_domain.clone().into_iter().map(|h|cp.evaluate(h)).collect();
let cp_merkle_tree=MerkleTree::new(&cp_evaluations);
//sending merkle root to channel i.e. verifier
channel.send(cp_merkle_tree.inner.root().unwrap().to_vec());
// now we have made a composite polynomial if p_0 , p_1 and p_2 are all polynomial constraints then cp is also a polynomial constraint 
//we can prove that using F.R.I. (Fast Reedman-Impaliacio) protocol
// we can prove composite polynomial is close to a polynomial by showing it's distance from the polynomial is small
// FRI is a folding scheme which convert the polynomial into even and odd polynomial of smaller degree
//prover tries to convience the verifier that the FRI folding is true by commitiming it to merkle tree
// Part 2: FRI
 // 1. Receive random element `beta` from verifier.
        // 2. Apply the FRI folding scheme or FRI operator to the composition polynomial.
        // 3. Commit to the new polynomial obtained after applying FRI operator.
        // 4. Send the new commitment to the verifier.
        // 5. Repeat step 1-4, until the polynomial degree is less than accepted degree in terms of security. in this case repeat till degree is 0.
        // 6. Prover sends the result to the verifier.
       
    // F.R.I Operator or folding scheme:
    // from proving: function is close to a polynomial of degree < D
    // to proving: new function is close to a new polynomial of degree < D/2, where new function has half the domain size of old polynomial.
    // Example: To prove: A function is close to a polynomial of a degree < 1024, with function domain size = 8192
    //          After applying the FRI operator we need to prove the new polynomial degree < 512 with new function domain size = 4096
    // split ot even and odd powers.
    // P_0(x) = g(x^2) + x h(x^2)
    // Get random element beta from verifier
    // P_1(y) = g(y) + beta * h(y)

    // For this example, repeat steps 1-4 till degree of polynomial < 1, when domain size is 8.

    // The new evaluation domain, will be half of the old evaluation domain.
    // and new evaluation domain is first half of the old evaluation domain squared.
    // eval domain: w, w.h, w.h2, .... w.h^8091
    // new eval domain: w^2, (w.h)^2, ... (w.h^4095)^2
    // square of the first half of the old eval domian, is equal to square of second half of old eval domain. This is a cyclic group property.

    // generate fri layers and commit to the fri layers.
    log::info!("generating fri layers and fri commitments");
    let (fri_polys, fri_domains, fri_layers, fri_merkle_trees) = fri_commit(
        cp,
        &eval_domain,
        &cp_evaluations,
        cp_merkle_tree.clone(),
        &mut channel,
    );

    assert!(fri_layers.len() == 11); // 11 fri layers, 8192 field elements in eval domain to 8 field elements in eval domain
    assert!(fri_layers[fri_layers.len() - 1].len() == 8); // last layer will have evaluations at 8 points
    assert!(fri_polys[fri_polys.len() - 1].degree() == 0); // and degree of last polynomial equal to zero.

    // Proof or provers work contains generating commitments and decommiting them, to convence the verifier over the integrity of the computation.
    // i)  Commitment ✅
    // ii) Decommitment -> query phase
    // Decommitment involves verifier sending random elements from evaluation domain to prover. and prover responding with decommitments to the evaluations, which involve sending merkle paths along with evaluations.
    //
    // there are a total of 12 commitments made and send over channel to verifier. The commitments made are: trace, composition polynomial, 10 fri layers.
    //
    // with each successful query and valid decommitment, verifiers confidence in the proof increases.

    log::info!("decommitting fri layers");
    let (num_of_queries, blow_up_factor, maximum_random_int) = (4, 8, 8192 - 16);
    decommit_fri(
        num_of_queries,
        blow_up_factor,
        8192 - 16,
        &f_evaluations,
        &f_merkle_tree,
        &fri_layers,
        &fri_merkle_trees,
        &mut channel,
    );

    log::info!(
        "proof generation complete, time taken: {}ms, proof size: {} bytes, compressed proof size: {} bytes",
        (Local::now() - start_time).num_milliseconds(),
        channel.proof_size(),
        channel.compressed_proof_size()
    );

    let compressed_proof = channel.compressed_proof;

    // verify the proof
    log::info!("verifying proof");
    let start = Local::now();
    verify_proof(
        num_of_queries,
        maximum_random_int,
        blow_up_factor,
        field,
        &fri_domains,
        &compressed_proof,
    );

    log::info!(
        "verification successful, time taken: {:?}µs",
        (Local::now() - start).num_microseconds().unwrap()
    );
} 















    










