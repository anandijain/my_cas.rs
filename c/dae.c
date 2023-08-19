#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

int residual(realtype t, N_Vector y, N_Vector ydot,
             N_Vector res, void *user_data)
{

    // Extract values from y vector
    realtype x = NV_Ith_S(y, 0);
    realtype y_val = NV_Ith_S(y, 1);

    // Extract derivatives from ydot vector
    realtype xdot = NV_Ith_S(ydot, 0);
    realtype ydot_val = NV_Ith_S(ydot, 1);

    // Parameters (e.g., gravity and pendulum length)
    realtype g = 9.81; // Gravity
    realtype L = 1.0;  // Pendulum length

    // Residuals for the DAE system
    NV_Ith_S(res, 0) = xdot - L * ydot_val;
    NV_Ith_S(res, 1) = ydot_val + L * xdot - g;

    // Algebraic constraint (y^2 = L^2 - x^2)
    NV_Ith_S(res, 2) = x * x + y_val * y_val - L * L;

    return 0;
}
int main()
{
    SUNContext sunctx = NULL; // Assuming NULL works as a default context. If not, proper initialization is required.

    void *mem = IDACreate(sunctx);

    N_Vector y = N_VNew_Serial(3, sunctx);    // Adjust size accordingly
    N_Vector ydot = N_VNew_Serial(3, sunctx); // Adjust size accordingly

    // Initial conditions
    NV_Ith_S(y, 0) = 0.0;
    NV_Ith_S(ydot, 0) = 0.0;

    realtype t = 0.0; // Define and initialize t
    IDAInit(mem, residual, t, y, ydot);
    IDASStolerances(mem, 1.0e-6, 1.0e-6);

    // Solve the DAE system, for example till t = 1.0
    realtype tout = 1.0;
    IDASolve(mem, tout, &t, y, ydot, IDA_NORMAL);

    // Cleanup
    N_VDestroy(y);
    N_VDestroy(ydot);
    IDAFree(&mem);

    return 0;
}