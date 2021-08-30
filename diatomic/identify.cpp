#include <iostream>
#include <iomanip>
// #include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


py::array_t<double> identify_by_largest_contribution(double **clevels, double **elevels,
                                                     int &nch, int &cshapex, int &eshapex)
{
    /* allocate the output buffer */
    int result_size = cshapex;
    if (eshapex > cshapex)
        result_size = eshapex;

    int ncols = 12+nch;  // number of columns in the output array
    py::array_t<double> result = py::array_t<double>(result_size*ncols);
    py::buffer_info result_buf = result.request();
    double *result_ptr = (double *) result_buf.ptr;

    int count = 0;
    int nrows = 0;
    for(int i = 0; i < eshapex; i++)
    {
        auto elevel = elevels[i];

        int eind = (int)elevel[0];
        int ev = (int) elevel[1];
        double ej = elevel[2];
        double eenrg = elevel[3];
        int ep = (int) elevel[4];
        double eunc = elevel[5];
        int em = (int) elevel[6];
        int es = (int) elevel[7];

        for(int j = 0; j < cshapex; j++)
        {
            auto clevel = clevels[j];

            int cind = (int) clevel[0];
            double cenrg = clevel[1];
            double cj = clevel[2];
            int cp = (int) clevel[3];
            int ci = (int) clevel[4];
            int cs = (int) clevel[5+nch];
            int cv = (int) clevel[6+nch];
            int cl = (int) clevel[7+nch];
            double co = clevel[8+nch];

            if (cv == ev && cj == ej && cp == ep && 
                cs == es && ci == ((em / 10) + 1))
            {
                result_ptr[count++] = cv;
                result_ptr[count++] = cj;
                result_ptr[count++] = co;
                result_ptr[count++] = co - cl;  // sigma 
                result_ptr[count++] = cl;
                result_ptr[count++] = cp;
                result_ptr[count++] = em;
                result_ptr[count++] = cenrg;
                result_ptr[count++] = eenrg;
                result_ptr[count++] = cenrg - eenrg;
                result_ptr[count++] = eunc;
                for (int k = 0; k < nch; k++)
                {
                    result_ptr[count++] = clevel[4+nch+k];
                }
                result_ptr[count++] = cs;
                nrows++;
                break;
            }
        }
    }

    /* reshape result */
    result.resize({nrows, ncols});

    return result;
}

py::array_t<double> identify_by_closest_energy(double **clevels, double **elevels,
                                               int &nch, int &cshapex, int &eshapex)
{
    /* allocate the output buffer */
    int result_size = cshapex;
    if (eshapex > cshapex)
        result_size = eshapex;

    int ncols = 12 + nch;  // number of columns in the output array
    py::array_t<double> result = py::array_t<double>(result_size*ncols);
    py::buffer_info result_buf = result.request();
    double *result_ptr = (double *) result_buf.ptr;

    int count = 0;
    int nrows = 0;
    int v_tmp = 0, lambda_tmp = 0, par_tmp = 0, marker_tmp = 0, state_tmp = 0;
    double j_tmp = 0.0, omega_tmp = 0.0, sigma_tmp = 0.0;
    double cenrg_tmp = 0.0, eenrg_tmp = 0.0;
    double diff_enrg_tmp = 0.0, unc_tmp = 0.0;
    double ccoefs_tmp [nch] = {0.0};

    // This algorithm works correctly if the calculated 
    // levels are sorted by energy in ascending order

    for(int i = 0; i < eshapex; i++)
    {
        auto elevel = elevels[i];

        int eind = (int)elevel[0];
        int ev = (int) elevel[1];
        double ej = elevel[2];
        double eenrg = elevel[3];
        int ep = (int) elevel[4];
        double eunc = elevel[5];
        int em = (int) elevel[6];
        int es = (int) elevel[7];

        double min_diff_enrg = 1.0e10;
        bool is_diff_enrg_min = false;
        bool is_last_iter = false;

        // TODO: to account for used levels

        for(int j = 0; j < cshapex; j++)
        {
            auto clevel = clevels[j];

            int cind = (int) clevel[0];
            double cenrg = clevel[1];
            double cj = clevel[2];
            int cp = (int) clevel[3];
            int ci = (int) clevel[4];
            int cs = (int) clevel[5+nch];
            int cv = (int) clevel[6+nch];
            int cl = (int) clevel[7+nch];
            double co = clevel[8+nch];

            double diff_enrg = std::abs(cenrg - eenrg);

            if (cj == ej && cp == ep && ci == ((em / 10) + 1))
            {
                if (diff_enrg < min_diff_enrg)
                {
                    is_diff_enrg_min = false;
                    min_diff_enrg = diff_enrg;

                    v_tmp = ev;
                    j_tmp = cj;
                    omega_tmp = co;
                    sigma_tmp = co - cl;  // sigma 
                    lambda_tmp = cl;
                    par_tmp = cp;
                    marker_tmp = em;
                    cenrg_tmp = cenrg;
                    eenrg_tmp = eenrg;
                    diff_enrg_tmp = cenrg - eenrg;
                    unc_tmp = eunc;
                    for (int k = 0; k < nch; k++)
                    {
                        ccoefs_tmp[k] = clevel[4+nch+k];
                    }
                    state_tmp = es;
                }
                else
                {
                    is_diff_enrg_min = true;
                    break;
                }
            }
            // when the correct identified level is the first met
            if (j == cshapex - 1)
            {
                is_last_iter = true;
            }
        }

        if (is_diff_enrg_min == true || is_last_iter == true)
        {
            result_ptr[count++] = v_tmp;
            result_ptr[count++] = j_tmp;
            result_ptr[count++] = omega_tmp;
            result_ptr[count++] = omega_tmp - lambda_tmp;  // sigma 
            result_ptr[count++] = lambda_tmp;
            result_ptr[count++] = par_tmp;
            result_ptr[count++] = marker_tmp;
            result_ptr[count++] = cenrg_tmp;
            result_ptr[count++] = eenrg_tmp;
            result_ptr[count++] = diff_enrg_tmp;
            result_ptr[count++] = unc_tmp;
            for (int k = 0; k < nch; k++)
            {
                result_ptr[count++] = ccoefs_tmp[k];
            }
            result_ptr[count++] = state_tmp;
            nrows++;
        }
    }

    /* reshape result */
    result.resize({nrows, ncols});

    return result;
}

py::array_t<double> identify_levels(py::array_t<double> clevels_inp, 
                                    py::array_t<double> elevels_inp,
                                    int &nch, int &identify_option)
{
    /* read input arrays buffer info */
    py::buffer_info clevels_buf = clevels_inp.request();
    py::buffer_info elevels_buf = elevels_inp.request();

    double *clevesls_ptr = (double *) clevels_buf.ptr;
    double *elevels_ptr = (double *) elevels_buf.ptr;

    int cshapex = clevels_buf.shape[0];
    int cshapey = clevels_buf.shape[1];
    int eshapex = elevels_buf.shape[0];
    int eshapey = elevels_buf.shape[1];
    int cl_shapex = cshapex * eshapex;
    int cl_shapey = 2*cshapey + 2;

    /* allocate 2D arrays and fill them */
    double **clevels = new double*[cshapex];
    for(int i = 0; i < cshapex; i++)
        clevels[i] = new double[cshapey];

    for (int idx = 0; idx < cshapex; idx++)
        for (int idy = 0; idy < cshapey; idy++)
            clevels[idx][idy] = clevesls_ptr[idx*cshapey+idy];

    double **elevels = new double*[eshapex];
    for(int i = 0; i < eshapex; i++)
        elevels[i] = new double[eshapey];

    for (int idx = 0; idx < eshapex; idx++)
        for (int idy = 0; idy < eshapey; idy++)
            elevels[idx][idy] = elevels_ptr[idx*eshapey+idy];


    if (identify_option == 0)
    {
        return identify_by_closest_energy(clevels, elevels, nch, cshapex, eshapex);
    }
    else //elseif (identify_option == 1)
    {
        return identify_by_largest_contribution(clevels, elevels, nch, cshapex, eshapex);
    }

    /* delete allocated memory */
    for(int i = 0; i < cshapex; i++)
        delete [] clevels[i];
    delete [] clevels;

    for(int i = 0; i < eshapex; i++)
        delete [] elevels[i];
    delete [] elevels;
}

PYBIND11_MODULE(identify, m)
{
    m.def("identify_levels", 
          &identify_levels, 
          "A function which identifies levels");
}
