#include <iostream>
#include <iomanip>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<double> calculate_wavenumbers_list(py::array_t<double> ulevels_inp, 
											   py::array_t<double> llevels_inp)
{
	/* read input arrays buffer info */
	py::buffer_info ulevels_buf = ulevels_inp.request();
	py::buffer_info llevels_buf = llevels_inp.request();

	double *ulevesls_ptr = (double *) ulevels_buf.ptr;
	double *llevels_ptr = (double *) llevels_buf.ptr;

	int ushapex = ulevels_buf.shape[0];
	int ushapey = ulevels_buf.shape[1];
	int lshapex = llevels_buf.shape[0];
	int lshapey = llevels_buf.shape[1];
	int cw_shapex = ushapex * lshapex;
	int cw_shapey = 2*ushapey + 2;

	/* allocate the output buffer */
	int result_size = cw_shapex;

	int ncols = 18;  // number of columns in the output array
	py::array_t<double> result = py::array_t<double>(result_size*ncols);
	py::buffer_info result_buf = result.request();

	double *result_ptr = (double *) result_buf.ptr;

	/* allocate 2D arrays and fill them */
	double **ulevels = new double*[ushapex];
	for(int i = 0; i < ushapex; i++)
		ulevels[i] = new double[ushapey];

	for (int idx = 0; idx < ushapex; idx++)
		for (int idy = 0; idy < ushapey; idy++)
			ulevels[idx][idy] = ulevesls_ptr[idx*ushapey+idy];

	double **llevels = new double*[lshapex];
	for(int i = 0; i < lshapex; i++)
		llevels[i] = new double[lshapey];

	for (int idx = 0; idx < lshapex; idx++)
		for (int idy = 0; idy < lshapey; idy++)
			llevels[idx][idy] = llevels_ptr[idx*lshapey+idy];

	/* obtain the list of calculated wavenumbers */
	int wcount = 0;
	int nrows = 0;
	for(int i = 0; i < ushapex; i++)
	{
		auto ulevel = ulevels[i];

		int uind = (int)ulevel[0];
		int uv = (int)ulevel[1];
		double uj = ulevel[2];
		double uenrg = ulevel[3];      
		int up = (int)ulevel[4];
		int us = (int)ulevel[5];
		int ulam = ulevel[7];
		int uomg = (int)ulevel[8];

		for(int j = 0; j < lshapex; j++)
		{
			auto llevel = llevels[j];

			int lind = (int)llevel[0];
			int lv = (int)llevel[1];
			double lj = llevel[2];
			double lenrg = llevel[3];
			int lp = (int)llevel[4];
			int ls = (int)llevel[5];
			int llam = llevel[7];
			int lomg = (int)llevel[8];

			// -1 means the transition is not allowed by the strict 
			// selection rules for J and symmetry
			int branch = -1;

			// selection rule for Pee lines
			if (uj == lj-1 && up == 1 && lp == 1){
				branch = 0;
			}
			// selection rule for Pff lines
			if (uj == lj-1 && up == 0 && lp == 0){
				branch = 1;
			}
			// selection rule for Qef lines
			if (uj == lj && up == 1 && lp == 0){
				branch = 2;
			}
			// selection rule for Qfe lines
			if (uj == lj && up == 0 && lp == 1){
				branch = 3;
			}
			// selection rule for Ree lines
			if (uj == lj+1 && up == 1 && lp == 1){
				branch = 4;
			}
			// selection rule for Rff lines
			if (uj == lj+1 && up == 0 && lp == 0){
				branch = 5;
			}

			double cwaven = uenrg - lenrg;

			if(cwaven > 0.0)
			{
				result_ptr[wcount++] = uind;
				result_ptr[wcount++] = uv;
				result_ptr[wcount++] = uj;
				result_ptr[wcount++] = up;
				result_ptr[wcount++] = us;
				result_ptr[wcount++] = uenrg;
				result_ptr[wcount++] = ulam;
				result_ptr[wcount++] = uomg;
				result_ptr[wcount++] = lind;
				result_ptr[wcount++] = lv;
				result_ptr[wcount++] = lj;
				result_ptr[wcount++] = lp;
				result_ptr[wcount++] = ls;
				result_ptr[wcount++] = lenrg;
				result_ptr[wcount++] = llam;
				result_ptr[wcount++] = lomg;
				result_ptr[wcount++] = cwaven;
				result_ptr[wcount++] = branch;
				nrows++;
			}
		}
	}

	/* reshape result */
	result.resize({nrows, ncols});

	/* delete allocated memory */
	for(int i = 0; i < ushapex; i++)
		delete [] ulevels[i];
	delete [] ulevels;

	for(int i = 0; i < lshapex; i++)
		delete [] llevels[i];
	delete [] llevels;

	return result;
}

PYBIND11_MODULE(wavenumbers, m)
{
	m.def("calculate_wavenumbers_list", 
		  &calculate_wavenumbers_list, 
		  "A function which calculates the final wavenumbers list");
}
