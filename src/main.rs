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
log::info!("commiting to LDE");
let f_merkle_tree=MerkleTree::new(&f_evaluations);
let mut channel = Channel::new();
channel.send(f_merkle_tree.inner.root().unwrap().to_vec());
log::info!("generating constraint polynomials");
log::debug!("generating p_0");
let p_0=(trace_polynomial.clone())-Polynomial::new_from_coefficients(vec![FieldElement(1,field)])/Polynomial::new_from_coefficients(vec![-g.pow(0),FieldElement(1,field)]);
let p_1=(trace_polynomial.clone())-Polynomial::new_from_coefficients(vec![FieldElement(2338775057,field)])/Polynomial::new_from_coefficients(vec![-g.pow(1022),FieldElement(1,field)]);
let mut  x_1024=Polynomial::new_from_coefficients(vec![FieldElement::new(0,field);1025]);
x_1024.coefficients[0]=-FieldElement::new(1,field);
x_1024.coefficients[1024]=FieldElement::new(1,field);
let p_2_den = x_1024/Polynomial::new_from_coefficients(vec![-g.pow(1021),FieldElement(1,field)])* Polynomial::new_from_coefficients(vec![-g.pow(1022),FieldElement(1,field)])*Polynomial::new_from_coefficients(vec![-g.pow(1023),FieldElement(1,field)]);
let trace_poly_g2=trace_polynomial.clone().compose(g.pow(2));
let trace_poly_g= trace_polynomial.clone().compose(g);
let square_trace_poly=trace_polynomial.clone()*trace_polynomial.clone();
let p2_num=trace_poly_g2-trace_poly_g-square_trace_poly;
let p_2= p2_num/p_2_den ;
log::info!("receive alpha_0,aplha_1,alpha_2");
let alpha_0=channel.receive_random_field_element(field);
let alpha_1=channel.receive_random_field_element(field);
let alpha_2=channel.receive_random_field_element(field);
let cp = Polynomial::new_from_coefficients(vec![alpha_0])*p_0 + Polynomial::new_from_coefficients(vec![alpha_1])*p_1 + Polynomial::new_from_coefficients(vec![alpha_2])*p_2;














    









}
