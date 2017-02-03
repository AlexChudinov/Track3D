#include <iostream>
#include <libIntegrators.h>
#include <fstream>

/*Simple harmonic oscillator:*/
struct diff_x2
{
	double a;
	double b;
	static void diff(const void* pObj, 
		const double * state, 
		double * dxdt, 
		const double * /*time*/)
	{
		const diff_x2* ptr = reinterpret_cast<const diff_x2*>(pObj);
		//dx/dt = -a*y
		//dy/dt = b*x
		//d2y/dt2 = -b*a*y => y = cos(sqrt(b*a)*y + phi0)
		dxdt[0] = - ptr->a * state[1];
		dxdt[1] = ptr->b * state[0];
	}
};


int main()
{
	diff_x2 Obj = { 1,1 }; //Set model parameters

	//Create integrator interface.
	void* pObj = create_integrator_interface
	(
		2,								//Set state size
		RUNGE_KUTTA4,			//Choose integrator type
		reinterpret_cast<void*>(&Obj),	//Put your model and parameters
		diff_x2::diff					//Put pointer to a function which calculates derivatives
	);

	if (!pObj) //Check if everything is ok
	{
		std::cout << "Failed to create integrator interface.\n";
		return 0;
	}
	else
	{
		std::cout << "Integrator interface was created.\n";
	}

	std::ofstream out;
	out.open("test.txt"); //Save result into the file

	double state[2] = { 0, 1 }; //Initialize state
	double dt = 0.5;			//Time-step
	for (double time = 0.0; time < 10; time += dt)
	{
		std::cout << state[0] << "\t" << state[1] << std::endl;
		out << time << "\t" << state[0] << "\t" << state[1] << std::endl;
		do_integrator_step(pObj, state, &time, &dt); //Proceed one integrator time-step
	}
	out.close();

	delete_integrator_interface(pObj); //Free integrator
	return 0;
}

