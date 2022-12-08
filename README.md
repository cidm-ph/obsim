# obsim

> **Outbreak simulation based on a branching process**

Simulate outbreaks of a communicable disease.
This is both a simulation framework and an implementation that uses very
simplistic modelling assumptions.

## How to use

This is a Rust library, but comes with some [examples](./examples) that you can modify.
See [Installing Rust](https://www.rust-lang.org/tools/install) to get set up.

With the repository checked out, you can run these examples with e.g.
`cargo run --example simple` to get an annotated FASTA file with the simulation result.

## Purpose

The currently implemented models are too simplistic to capture many real
features of microbial evolution. They are intended for simple approximations
that apply on short timescales where only a few mutations are expected to occur
over the duration of outbreaks.

The framework can be extended with more sophisticated models using the same
simple interfaces to drive simulations.

## Citation

This framework was created for the following study

> Suster CJE, Arnott A, Blackwell G, Gall M, Draper J, Martinez E, Drew AP, Rockett RJ, Chen SC-A, Kok J, Dwyer DE and Sintchenko V (2022)
> Guiding the design of SARS-CoV-2 genomic surveillance by estimating the resolution of outbreak detection.
> Front. Public Health 10:1004201. doi: [10.3389/fpubh.2022.1004201](https://doi.org/10.3389/fpubh.2022.1004201)

## Licence

Dual-licensed under [MIT](LICENSE-MIT) or [Apache 2.0](LICENSE-APACHE).

© 2022 Western Sydney Local Health District, NSW Health
