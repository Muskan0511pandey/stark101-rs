use crate::channel::Channel;
use crate::field::{Field, FieldElement};
use crate::merkle_tree::MerkleTree;
///
/// panics for invalid proof.
pub fn verify_proof(
    num_of_queries: usize,
    maximum_random_int: u64,
    blow_up_factor: usize,
    field: Field,
    fri_domains: &[Vec<FieldElement>],
    compressed_proof: &[Vec<u8>],
) {
    let mut channel = Channel::new();

    // merkle root of trace polynomial LDE
    let f_merkle_root = compressed_proof[0].clone();
    channel.send(f_merkle_root.clone());

    // get random values for calculating composition polynomial
    let _alpha_0 = channel.receive_random_field_element(field);
    let _alpha_1 = channel.receive_random_field_element(field);
    let _alpha_2 = channel.receive_random_field_element(field);

    let cp_merkle_root = compressed_proof[1].clone();
    channel.send(cp_merkle_root.clone());

    // commit to fri
    let mut fri_merkle_roots = Vec::new();
    let mut betas = Vec::new();
    for i in 0..10 {
        let beta = channel.receive_random_field_element(field);
        betas.push(beta);
        let fri_root = compressed_proof[2 + i].clone();
        channel.send(fri_root.clone());
        fri_merkle_roots.push(fri_root);
    }

    let last_layer_free_term = compressed_proof[12].clone();
    channel.send(last_layer_free_term.clone());

    let mut base_idx = 13;

    // decommit fri
    for i in 0..num_of_queries {
        let idx = channel.receive_random_int(0, maximum_random_int, true) as usize;
        // verifiy the queries too.
        verify_queries(
            base_idx + i,
            idx,
            blow_up_factor,
            field,
            &fri_merkle_roots,
            &fri_domains,
            compressed_proof,
            &betas,
            &mut channel,
        );
        base_idx += 46;
    }
}

pub fn verify_fri_layers(
    base_idx: usize,
    idx: usize,
    field: Field,
    fri_merkle_roots: &[Vec<u8>],
    fri_domains: &[Vec<FieldElement>],
    compressed_proof: &[Vec<u8>],
    betas: &[FieldElement],
    channel: &mut Channel,
) {
    let lengths = vec![8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16];
    for i in 0..10 {
        // verify the fri layers
        let length = lengths[i];
        let elem_idx = idx % length;
        let elem = compressed_proof[base_idx + 4 * i].clone();
        channel.send(elem.clone());
        let elem_proof = compressed_proof[base_idx + 4 * i + 1].clone();
        channel.send(elem_proof.clone());
        let merkle_root = if i == 0 {
            compressed_proof[1].clone()
        } else {
            fri_merkle_roots[i - 1].clone()
        };
        if i != 0 {
            // checking fri polynomials consistency.
            let prev_elem =
                FieldElement::from_bytes(&compressed_proof[base_idx + 4 * (i - 1)].clone());
            let prev_sibling =
                FieldElement::from_bytes(&compressed_proof[base_idx + 4 * (i - 1) + 2].clone());
            let two = FieldElement::new(2, field);
            let computed_elem = (prev_elem + prev_sibling) / two
                + (betas[i - 1] * (prev_elem - prev_sibling)
                    / (two * fri_domains[i - 1][idx % lengths[i - 1]]));
            assert!(computed_elem.0 == FieldElement::from_bytes(&elem).0);
        }
        assert!(MerkleTree::validate(
            merkle_root.clone(),
            elem_proof,
            elem_idx,
            elem.clone(),
            length,
        ));
        let sibling_idx = (idx + length / 2) % length;
        let sibling = compressed_proof[base_idx + 4 * i + 2].clone();
        channel.send(sibling.clone());
        let sibling_proof = compressed_proof[base_idx + 4 * i + 3].clone();
        channel.send(sibling_proof.clone());

        assert!(MerkleTree::validate(
            merkle_root,
            sibling_proof,
            sibling_idx,
            sibling.clone(),
            length,
        ));
    }
    let last_elem = compressed_proof[base_idx + 40].clone();
    channel.send(last_elem);
}

pub fn verify_queries(
    base_idx: usize,
    idx: usize,
    blow_up_factor: usize,
    field: Field,
    fri_merkle_roots: &[Vec<u8>],
    fri_domains: &[Vec<FieldElement>],
    compressed_proof: &[Vec<u8>],
    betas: &[FieldElement],
    channel: &mut Channel,
) {
    let len = 8192;
    let f_merkle_root = compressed_proof[0].clone();
    // verifiy the queries too.
    let f_x = compressed_proof[base_idx].clone();
    channel.send(f_x.clone());
    let f_x_auth = compressed_proof[base_idx + 1].clone();
    channel.send(f_x_auth.clone());

    assert!(MerkleTree::validate(
        f_merkle_root.clone(),
        f_x_auth,
        idx,
        f_x,
        len
    ));
    let f_gx = compressed_proof[base_idx + 2].clone();
    channel.send(f_gx.clone());
    let f_gx_auth = compressed_proof[base_idx + 3].clone();
    channel.send(f_gx_auth.clone());
    assert!(MerkleTree::validate(
        f_merkle_root.clone(),
        f_gx_auth,
        idx + blow_up_factor,
        f_gx,
        len
    ));
    let f_g2x = compressed_proof[base_idx + 4].clone();
    channel.send(f_g2x.clone());
    let f_g2x_auth = compressed_proof[base_idx + 5].clone();
    channel.send(f_g2x_auth.clone());
    assert!(MerkleTree::validate(
        f_merkle_root.clone(),
        f_g2x_auth,
        idx + 2 * blow_up_factor,
        f_g2x,
        len
    ));
    // can we actually able to relate cp and f_g2x,f_gx, f_x?
    // it looks like we need to generate p1, p2 and p3 for relating  cp and f_g2x, f_gx, f_x -> so, verifier could not do that.
    verify_fri_layers(
        base_idx + 6,
        idx,
        field,
        fri_merkle_roots,
        fri_domains,
        compressed_proof,
        betas,
        channel,
    );
}
