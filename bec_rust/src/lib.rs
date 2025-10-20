use fftw::{
    types::{
        Sign,
        c64
    },
    plan,
};

use num::{Complex};
use lapack::*;

// Type Definitions for pure convenience (and for making) public maybe.
type Direction = Sign;
type Plan = plan::C2CPlan64;

// Somee usevul Constants
pub const PI: f64 = std::f64::consts::PI;
pub const FRAC_ROOT_TWO_PI: f64 = 0.398942280401432677939946059934381868_f64;
pub const E: f64 = std::f64::consts::E;

const I: c64 = Complex::I;

pub mod physics {

    use num::pow;

    fn potential(location: f64, &wave_number: &f64, trap: bool, lattice: bool) -> f64 {
        match (trap, lattice) {
            (true, false) => pow(location, 2) + 0.5,
            (false, true) => {
                let sinx = f64::sin( wave_number * location );
                let pot = 0.5 * pow(sinx, 2) * pow (location, 2);
                pot
            },
            (true, true) => {
                let sinx = f64::sin( wave_number * location );
                let lattice = 0.5 * pow(sinx, 2) * pow (location, 2);
                let pot = 0.5 * pow(sinx, 2) * pow (location, 2) + 0.5 * pow(sinx, 2) * pow (location, 2);
                pot
            },
            _ => 0.
        }
   }
}

/// Linear Algebra used.
/// at the moment, the eigenvalues and eigenvectors of a tridiagonal Matrix are used,
/// so we will use LAPACK's dsyev
pub mod linalg {

    pub enum Jobz {
        EigenValuesOnly,
        WithEigenvectors
    }
    fn get_jobz(selection: Jobz) -> char {
        match selection {
            Jobz::EigenValuesOnly => 'N',
            Jobz::WithEigenvectors => 'V',
        }
    }

    pub enum Uplo {
        UpperTriangle,
        LowerTriangle
    }
    fn get_uplo(upperOrLower: Uplo) -> char {
        match upperOrLower {
            Uplo::UpperTriangle => 'U',
            Uplo::LowerTriangle => 'L',
        }
    }

    /// #Config
    ///
    /// Configuration struct for LAPACK Functions in General
    /// Parameters:
    /// * `jobz`: This has been put into an enum for readability
    /// * `n`: Matrix order. usize. We only need one rank, since we want a symmetric matrix
    pub struct EigenConfig {
        jobz: char,
        uplo: char,
        n: i32,
        system_width: f64,
    }
    impl EigenConfig {
        pub fn init(job_size: Jobz, upper_lower: Uplo, step_number: usize, system_width: f64 ) -> EigenConfig {
            let jobz = get_jobz(job_size);
            let uplo = get_uplo(upper_lower);
            let n = step_number as i32;
            let system_width = system_width;

            EigenConfig { jobz, uplo, n, system_width }
        }
    }
}
