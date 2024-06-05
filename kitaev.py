from mpmath import re, jtheta, exp, acos, sqrt
from mpmath import j as i
from qat.lang import qrout, RY, X


def JacobiTheta(mean, var):
    """Calculates an ad-hoc Jacobi theta function as an mpmath floating point number

    Parameters
    ----------
    mean : float
        Gaussian state mean parameter
    var : float
        Gaussian state variance parameter

    Returns
    -------
    mpmath.ctx_mp_python.mpf
        Jacobi theta function value
    """

    return re(jtheta(3, mean/(i*var**2), exp(-1/var**2))) 

def RotationAngle(mean, var):
    """Calculates Kitaev's rotation angle for Gaussian state preparation

    Parameters
    ----------
    mean : float
        Gaussian function mean parameter
    var : float
        Gaussian function variance parameter

    Returns
    -------
    float
        Quantum gate rotation angle
    """

    cos_alpha_sqrd = JacobiTheta(mean/2, var/2) / JacobiTheta(mean, var)
    alpha = acos(sqrt(cos_alpha_sqrd))

    return float(alpha)


@qrout
def GaussianState(mean, var, qbits):
    """Create Kitaev's quantum gate for preparing a discretized single-variate Gaussian state

    Parameters
    ----------
    mean : float
        Gaussian function mean parameter
    var : float
        Gaussian function variance parameter
    qbits : int
        Number of qbits allocated for discretization

    Returns
    -------
    qat.lang.AQASM.routines.QRoutine
        Kitaev's Gaussian state preparation gate
    """

    # Rotate rightmost qubit
    RY(2*RotationAngle(mean,var))(qbits-1) 
    # The actual state preparation subroutine is defined recursively
    if qbits>1:
        # Switch rightmost bit logical values
        X(qbits-1) 
        # Complete one half of the state
        GaussianState(mean/2, var/2, qbits-1).ctrl()(qbits-1,range(qbits-1)) 
        # Unswitch rightmost bit logical values
        X(qbits-1)
        # Complete the other half of the state
        GaussianState((mean-1)/2, var/2, qbits-1).ctrl()(qbits-1,range(qbits-1)) 