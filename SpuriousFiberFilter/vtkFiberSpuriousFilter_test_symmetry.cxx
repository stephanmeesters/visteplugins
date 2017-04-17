
/** Includes */

#include "vtkFiberSpuriousFilter.h"
#include <ctime>
#include <algorithm>
#include <omp.h>


namespace bmia {


vtkStandardNewMacro(vtkFiberSpuriousFilter);


//-----------------------------[ Constructor ]-----------------------------\\

vtkFiberSpuriousFilter::vtkFiberSpuriousFilter()
{
	// Set default options
    this->inputVolume = NULL;
    this->orientations =
        {
            0.,0.,1.,
        -0.20318274393026847,0.14762090442216314,0.9679487802288661,
        -0.42532540417602,0.3090169943749474,0.8506508083520399,
        -0.6095482317908053,0.4428627132664893,0.6575131712132815,
        -0.7236067977499789,0.5257311121191336,0.4472135954999579,
        0.07760890225389615,0.2388556408050596,0.9679487802288659,
        -0.1381966011250105,0.42532540417601994,0.8944271909999157,
        -0.36180339887498947,0.587785252292473,0.7236067977499789,
        -0.531939329536909,0.6817183540715489,0.5022953667054891,
        0.16245984811645317,0.5,0.8506508083520399,
        -0.052786404500042065,0.6881909602355867,0.7236067977499789,
        -0.2628655560595668,0.8090169943749473,0.5257311121191336,
        0.23282670676168846,0.7165669224151787,0.6575131712132815,
        0.02964396283142001,0.8641878268373419,0.5022953667054891,
        0.27639320225002106,0.8506508083520399,0.4472135954999579,
        -0.20318274393026847,-0.14762090442216314,0.9679487802288661,
        -0.42532540417602,-0.3090169943749474,0.8506508083520399,
        -0.6095482317908053,-0.4428627132664893,0.6575131712132815,
        -0.7236067977499789,-0.5257311121191336,0.4472135954999579,
        -0.4472135954999579,0.,0.8944271909999157,
        -0.6708203932499369,-0.16245984811645314,0.7236067977499789,
        -0.8127309757210738,-0.2952418088443262,0.5022953667054891,
        -0.6708203932499369,0.16245984811645314,0.7236067977499789,
        -0.85065080835204,0.,0.5257311121191336,
        -0.8127309757210738,0.2952418088443262,0.5022953667054891,
        0.07760890225389615,-0.2388556408050596,0.9679487802288659,
        0.16245984811645317,-0.5,0.8506508083520399,
        0.23282670676168846,-0.7165669224151787,0.6575131712132815,
        0.27639320225002106,-0.8506508083520399,0.4472135954999579,
        -0.1381966011250105,-0.42532540417601994,0.8944271909999157,
        -0.052786404500042065,-0.6881909602355867,0.7236067977499789,
        0.02964396283142001,-0.8641878268373419,0.5022953667054891,
        -0.36180339887498947,-0.587785252292473,0.7236067977499789,
        -0.2628655560595668,-0.8090169943749473,0.5257311121191336,
        -0.531939329536909,-0.6817183540715489,0.5022953667054891,
        0.2511476833527446,0.,0.9679487802288661,
        0.5257311121191336,0.,0.8506508083520399,
        0.7534430500582336,0.,0.6575131712132815,
        0.8944271909999159,0.,0.4472135954999579,
        0.36180339887498947,-0.2628655560595668,0.8944271909999157,
        0.6381966011250104,-0.2628655560595668,0.7236067977499789,
        0.8310519523121299,-0.23885564080505955,0.5022953667054891,
        0.44721359549995787,-0.5257311121191336,0.7236067977499789,
        0.6881909602355868,-0.5,0.5257311121191336,
        0.48397439011443294,-0.7165669224151788,0.5022953667054891,
        0.36180339887498947,0.2628655560595668,0.8944271909999157,
        0.44721359549995787,0.5257311121191336,0.7236067977499789,
        0.48397439011443294,0.7165669224151788,0.5022953667054891,
        0.6381966011250104,0.2628655560595668,0.7236067977499789,
        0.6881909602355868,0.5,0.5257311121191336,
        0.8310519523121299,0.23885564080505955,0.5022953667054891,
        0.7236067977499789,-0.5257311121191336,-0.4472135954999579,
        0.6095482317908053,-0.4428627132664893,-0.6575131712132815,
        0.42532540417602,-0.3090169943749474,-0.8506508083520399,
        0.20318274393026847,-0.14762090442216314,-0.9679487802288661,
        0.,0.,-1.,
        0.531939329536909,-0.6817183540715489,-0.5022953667054891,
        0.36180339887498947,-0.587785252292473,-0.7236067977499789,
        0.1381966011250105,-0.42532540417601994,-0.8944271909999157,
        -0.07760890225389615,-0.2388556408050596,-0.9679487802288659,
        0.2628655560595668,-0.8090169943749473,-0.5257311121191336,
        0.052786404500042065,-0.6881909602355867,-0.7236067977499789,
        -0.16245984811645317,-0.5,-0.8506508083520399,
        -0.02964396283142001,-0.8641878268373419,-0.5022953667054891,
        -0.23282670676168846,-0.7165669224151787,-0.6575131712132815,
        -0.27639320225002106,-0.8506508083520399,-0.4472135954999579,
        0.7236067977499789,0.5257311121191336,-0.4472135954999579,
        0.6095482317908053,0.4428627132664893,-0.6575131712132815,
        0.42532540417602,0.3090169943749474,-0.8506508083520399,
        0.20318274393026847,0.14762090442216314,-0.9679487802288661,
        0.8127309757210738,0.2952418088443262,-0.5022953667054891,
        0.6708203932499369,0.16245984811645314,-0.7236067977499789,
        0.4472135954999579,0.,-0.8944271909999157,
        0.85065080835204,0.,-0.5257311121191336,
        0.6708203932499369,-0.16245984811645314,-0.7236067977499789,
        0.8127309757210738,-0.2952418088443262,-0.5022953667054891,
        -0.27639320225002106,0.8506508083520399,-0.4472135954999579,
        -0.23282670676168846,0.7165669224151787,-0.6575131712132815,
        -0.16245984811645317,0.5,-0.8506508083520399,
        -0.07760890225389615,0.2388556408050596,-0.9679487802288659,
        -0.02964396283142001,0.8641878268373419,-0.5022953667054891,
        0.052786404500042065,0.6881909602355867,-0.7236067977499789,
        0.1381966011250105,0.42532540417601994,-0.8944271909999157,
        0.2628655560595668,0.8090169943749473,-0.5257311121191336,
        0.36180339887498947,0.587785252292473,-0.7236067977499789,
        0.531939329536909,0.6817183540715489,-0.5022953667054891,
        -0.8944271909999159,0.,-0.4472135954999579,
        -0.7534430500582336,0.,-0.6575131712132815,
        -0.5257311121191336,0.,-0.8506508083520399,
        -0.2511476833527446,0.,-0.9679487802288661,
        -0.8310519523121299,0.23885564080505955,-0.5022953667054891,
        -0.6381966011250104,0.2628655560595668,-0.7236067977499789,
        -0.36180339887498947,0.2628655560595668,-0.8944271909999157,
        -0.6881909602355868,0.5,-0.5257311121191336,
        -0.44721359549995787,0.5257311121191336,-0.7236067977499789,
        -0.48397439011443294,0.7165669224151788,-0.5022953667054891,
        -0.48397439011443294,-0.7165669224151788,-0.5022953667054891,
        -0.44721359549995787,-0.5257311121191336,-0.7236067977499789,
        -0.36180339887498947,-0.2628655560595668,-0.8944271909999157,
        -0.6881909602355868,-0.5,-0.5257311121191336,
        -0.6381966011250104,-0.2628655560595668,-0.7236067977499789,
        -0.8310519523121299,-0.23885564080505955,-0.5022953667054891,
        0.1552178045077923,0.9554225632202383,0.25114768335274457,
        -0.1381966011250105,0.9510565162951535,0.276393202250021,
        -0.4472135954999579,0.8506508083520399,0.276393202250021,
        -0.6871571340447014,0.6817183540715489,0.25114768335274457,
        0.,1.,0.,
        -0.3090169943749474,0.9510565162951535,0.,
        -0.5877852522924731,0.8090169943749473,0.,
        -0.1552178045077923,0.9554225632202383,-0.25114768335274457,
        -0.4360094506919568,0.8641878268373419,-0.25114768335274457,
        -0.8606959151435498,0.4428627132664894,0.2511476833527446,
        -0.9472135954999581,0.1624598481164532,0.2763932022500211,
        -0.9472135954999581,-0.1624598481164532,0.2763932022500211,
        -0.8606959151435498,-0.4428627132664894,0.2511476833527446,
        -0.9510565162951535,0.3090169943749474,0.,
        -0.9999999999999999,0.,0.,
        -0.9510565162951535,-0.3090169943749474,0.,
        -0.9566257939885021,0.1476209044221631,-0.25114768335274457,
        -0.9566257939885021,-0.1476209044221631,-0.25114768335274457,
        -0.6871571340447014,-0.6817183540715489,0.25114768335274457,
        -0.4472135954999579,-0.8506508083520399,0.276393202250021,
        -0.1381966011250105,-0.9510565162951535,0.276393202250021,
        0.1552178045077923,-0.9554225632202383,0.25114768335274457,
        -0.5877852522924731,-0.8090169943749473,0.,
        -0.3090169943749474,-0.9510565162951535,0.,
        0.,-1.,0.,
        -0.4360094506919568,-0.8641878268373419,-0.25114768335274457,
        -0.1552178045077923,-0.9554225632202383,-0.25114768335274457,
        0.4360094506919568,-0.8641878268373419,0.25114768335274457,
        0.6708203932499369,-0.6881909602355867,0.276393202250021,
        0.8618033988749894,-0.42532540417601994,0.276393202250021,
        0.9566257939885021,-0.1476209044221631,0.25114768335274457,
        0.5877852522924731,-0.8090169943749473,0.,
        0.8090169943749475,-0.587785252292473,0.,
        0.9510565162951535,-0.3090169943749474,0.,
        0.6871571340447014,-0.6817183540715489,-0.25114768335274457,
        0.8606959151435498,-0.4428627132664894,-0.2511476833527446,
        0.9566257939885021,0.1476209044221631,0.25114768335274457,
        0.8618033988749894,0.42532540417601994,0.276393202250021,
        0.6708203932499369,0.6881909602355867,0.276393202250021,
        0.4360094506919568,0.8641878268373419,0.25114768335274457,
        0.9510565162951535,0.3090169943749474,0.,
        0.8090169943749475,0.587785252292473,0.,
        0.5877852522924731,0.8090169943749473,0.,
        0.8606959151435498,0.4428627132664894,-0.2511476833527446,
        0.6871571340447014,0.6817183540715489,-0.25114768335274457,
        0.4472135954999579,-0.8506508083520399,-0.276393202250021,
        0.1381966011250105,-0.9510565162951535,-0.276393202250021,
        0.3090169943749474,-0.9510565162951535,0.,
        0.9472135954999581,0.1624598481164532,-0.2763932022500211,
        0.9472135954999581,-0.1624598481164532,-0.2763932022500211,
        1.,0.,0.,
        0.1381966011250105,0.9510565162951535,-0.276393202250021,
        0.4472135954999579,0.8506508083520399,-0.276393202250021,
        0.3090169943749474,0.9510565162951535,0.,
        -0.8618033988749894,0.42532540417601994,-0.276393202250021,
        -0.6708203932499369,0.6881909602355867,-0.276393202250021,
        -0.8090169943749475,0.587785252292473,0.,
        -0.6708203932499369,-0.6881909602355867,-0.276393202250021,
        -0.8618033988749894,-0.42532540417601994,-0.276393202250021,
        -0.8090169943749475,-0.587785252292473,0.
        };
}


//------------------------------[ Destructor ]-----------------------------\\

vtkFiberSpuriousFilter::~vtkFiberSpuriousFilter()
{

}

//-------------------------------[ SetInputVolume ]-------------------------------\\

void vtkFiberSpuriousFilter::SetInputVolume(vtkImageData * image)
{
	// Store the pointer
	this->inputVolume = image;
}

//-------------------------------[ SetParameters ]-------------------------------\\

void vtkFiberSpuriousFilter::SetParameters(ParameterSettings* ps)
{
    // Store the pointer
    this->ps = ps;
}

//-------------------------------[ Execute ]-------------------------------\\

double* Difference(double* vec, double* vec2)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = vec[0] - vec2[0];
    dp[1] = vec[1] - vec2[1];
    dp[2] = vec[2] - vec2[2];
    return dp;
}

double* Flip(double* vec)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = -vec[0];
    dp[1] = -vec[1];
    dp[2] = -vec[2];
    return dp;
}

double* HalvedDifference(double* vec, double* vec2)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = (vec[0] - vec2[0])/2.0;
    dp[1] = (vec[1] - vec2[1])/2.0;
    dp[2] = (vec[2] - vec2[2])/2.0;
    return dp;
}

double* Normalize(double* vec)
{
    double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    double* nvec = (double*) malloc(3*sizeof(double));
    nvec[0] = vec[0]/length;
    nvec[1] = vec[1]/length;
    nvec[2] = vec[2]/length;
    return nvec;
}

double* Cross(double* vec, double* vec2)
{
    double* c = (double*) malloc(3*sizeof(double));
    c[0] = vec[1]*vec2[2] - vec[2]*vec2[1];
    c[1] = vec[2]*vec2[0] - vec[0]*vec2[2];
    c[2] = vec[0]*vec2[1] - vec[1]*vec2[0];
    return c;
}

double Norm(double* vec)
{
    double output = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
    free(vec);
    return output;
}

void PrintVector3(double* vec)
{
    if(vec == NULL)
        return;
    printf("%f, %f, %f \n", vec[0], vec[1], vec[2]);
}

void PrintVector2(double* vec)
{
    if(vec == NULL)
        return;
    printf("{%f, %f} \n", vec[0], vec[1]);
}

void PrintMatrix3x3(double* mat)
{
    printf("{%f, %f, %f \n", mat[0], mat[1], mat[2]);
    printf(" %f, %f, %f \n", mat[3], mat[4], mat[5]);
    printf(" %f, %f, %f} \n", mat[6], mat[7], mat[8]);
}


double* EulerAngles(double* input)
{
    double x = input[0];
    double y = input[1];
    double z = input[2];
    double* output = (double*)malloc(2*sizeof(double));
    if((x - 1)*(x - 1) < 1e-10)
    {
        output[0] = 1.5708;
        output[1] = 0;
    }
    else if((x + 1)*(x + 1) < 1e-10)
    {
        output[0] = -1.5708;
        output[1] = 0;
    }
    else
    {
        output[0] = acos(sqrt(y*y + z*z)*(z>=0?1:-1))*(x>0?1:-1);
        output[1] = -asin(y/sqrt(y*y + z*z)*(z>=0?1:-1));
    }
    return output;
}

double* R(double* input)
{
    double beta = input[0];
    double gamma = input[1];
    double* output = (double*)malloc(9*sizeof(double));
    output[0] = cos(beta);
    output[1] = 0;
    output[2] = sin(beta);
    output[3] = sin(beta)*sin(gamma);
    output[4] = cos(gamma);
    output[5] = -cos(beta)*sin(gamma);
    output[6] = -cos(gamma)*sin(beta);
    output[7] = sin(gamma);
    output[8] = cos(beta)*cos(gamma);

    free(input);

    return output;
}

double* Transpose3x3(double *src)
{
    double* dst = (double*)malloc(sizeof(double)*9);
    for(int n = 0; n<9; n++)
    {
        int i = n/3;
        int j = n%3;
        dst[n] = src[3*j + i];
    }

    free(src);

    return dst;
}

double* Multiply(double* mat, double* vec)
{
    double* dst = (double*)malloc(sizeof(double)*3);
    dst[0] = mat[0]*vec[0] + mat[1]*vec[1]+mat[2]*vec[2];
    dst[1] = mat[3]*vec[0] + mat[4]*vec[1]+mat[5]*vec[2];
    dst[2] = mat[6]*vec[0] + mat[7]*vec[1]+mat[8]*vec[2];

    free(mat);

    return dst;
}

double* Subtract3(double* vec, double* vec2)
{
    double* dst = (double*)malloc(sizeof(double)*3);
    dst[0] = vec[0]-vec2[0];
    dst[1] = vec[1]-vec2[1];
    dst[2] = vec[2]-vec2[2];
    return dst;
}

double cot(double i) { return(1 / tan(i)); }

double kernel(double x, double y, double z, double b, double g, double D33, double D44, double t)
{
    return
    10.026513098524001*sqrt(D33*D44)*t*sqrt(D33*t)
    *
    (abs(b) < PI
    ?
    1.0/(100.53096491487338*D33*D44*t*t)*
        exp((-0.25*sqrt((1.*pow(-0.25*b*z +
                   (x*cos(0.5*b))/(1. - 0.041666666666666664*b*b),2))/(D33*D44)\
               + pow(b*b/D44 +
                pow(0.5*b*x + (0.5*z*cos(0.5*b))/
                    (1. - 0.041666666666666664*pow(b,2)),2)/D33,2)))/t)
    :
    1./(100.53096491487338*D33*D44*t*t)*
        exp((-0.25*sqrt((1.*pow(-0.25*b*z + 0.5*b*x*cot(0.5*b),2))/(D33*D44) +
              pow((b*b)/D44 + pow(0.5*b*x + 0.25*b*z*cot(0.5*b),2)/D33,2)))/t)
        )
    *
    (abs(g) < 0.3141592653589793
    ?
       1./(100.53096491487338*D33*D44*t*t)*
        exp((-0.25*sqrt((1.*pow(-0.25*g*z -
                   (1.*y*cos(0.5*g))/(1. - 0.041666666666666664*pow(g,2)),2))/
               (D33*D44) + pow((g*g)/D44 +
                pow(-0.5*g*y + (0.5*z*cos(0.5*g))/
                    (1. - 0.041666666666666664*pow(g,2)),2)/D33,2)))/t)
    :
       1./(100.53096491487338*D33*D44*pow(t,2))*
        exp((-0.25*sqrt((1.*pow(-0.25*g*z - 0.5*g*y*cot(0.5*g),2))/(D33*D44) +
              pow(pow(g,2)/D44 + pow(-0.5*g*y + 0.25*g*z*cot(0.5*g),2)/D33,2)))/t
          )
    );

}



double vtkFiberSpuriousFilter::k2(double* x, double* y, double* r, double* v)
{
    //printf("x:{%f,%f,%f}, y:{%f,%f,%f}, r:{%f,%f,%f}, v:{%f,%f,%f} \n", x[0], x[1], x[2], y[0], y[1], y[2], r[0], r[1], r[2], v[0], v[1], v[2]);

    //double* arg11 = R(EulerAngles(v));
    //PrintMatrix3x3(arg11);

    double* a = Subtract3(x,y);
    double* arg1 = Multiply(Transpose3x3(R(EulerAngles(v))),a);
    //PrintVector3(arg1);
//    PrintVector3(Subtract3(x,y));
//    PrintMatrix3x3(Transpose3x3(R(EulerAngles(v))));
//    PrintMatrix3x3(R(EulerAngles(v)));
//    PrintVector2(EulerAngles(v));

    //double* arg2 = EulerAngles(Multiply(Transpose3x3(R(EulerAngles(v))),r));
    double* arg2p = Multiply(Transpose3x3(R(EulerAngles(v))),r);
    double* arg2 = EulerAngles(arg2p);


//    PrintMatrix3x3(Transpose3x3(R(EulerAngles(v))));
//    PrintVector3(r);
//    PrintVector3(Multiply(Transpose3x3(R(EulerAngles(v))),r));
//    PrintVector2(arg2);

    //printf("{%f,%f,%f,%f,%f}\n",arg1[0],arg1[1],arg1[2],arg2[0],arg2[1]);





    double kernelval = kernel(arg1[0],arg1[1],arg1[2],arg2[0],arg2[1], this->ps->D33, this->ps->D44, this->ps->t);

    free(arg1);
    free(a);
    free(arg2);
    free(arg2p);

    //sprintf("KERNEL VALUE: %f \n",kernelval);
    return kernelval;


    //double aap[3];
   //Transpose3x3(arg1,aap,3,3);


    //PrintMatrix3x3(aap);

//    double* arg1 = MatrixMult(Transpose(R(EulerAngles(v))),Subtract(x,y));
//    double* arg2 = EulerAngles(MatrixMult(Transpose(R(EulerAngles(v)))));
//    return kernel(arg1[0],arg1[1],arg1[2],arg2[0],arg2[1]);
}


void vtkFiberSpuriousFilter::Execute()
{
    //
    //      CHECK POLYDATA INPUT
    //

	// Get the input
	vtkPolyData * input = this->GetInput();
	if(!input)
	{
		vtkErrorMacro(<< "Input has not been set.");
		return;
	}

	// Check if the input contains point data
	vtkPointData * inputPD = input->GetPointData();
	if (!inputPD)
	{
		vtkErrorMacro(<< "Input does not have point data.");
		return;
	}

	// Get the points of the input
	vtkPoints * inputPoints = input->GetPoints();
	if (!inputPoints)
	{
		vtkErrorMacro(<< "Input does not have points.");
		return;
	}

	// Get the lines array of the input
	vtkCellArray * inputLines = input->GetLines();
	if (!inputLines)
	{
		vtkErrorMacro(<< "Input does not have lines.");
		return;
	}

	// Get the output
	vtkPolyData * output = this->GetOutput();
	if (!output)
	{
		vtkErrorMacro(<< "Output has not been set.");
		return;
	}

	// Check if the output contains point data
	vtkPointData * outputPD = output->GetPointData();
	if (!outputPD)
	{
		vtkErrorMacro(<< "Output does not have point data.");
		return;
	}

    //
    //      PREPARE OUTPUT POLYDATA
    //

    // Prepare scalars list
    QList<vtkDoubleArray*> outputScalarsList;
    int numberOfScalarTypes = inputPD->GetNumberOfArrays();
    for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = vtkDoubleArray::New();
        outputScalars->SetName(inputPD->GetArray(i)->GetName());
        outputScalarsList.append(outputScalars);
    }

    // Add new scalar list for SM
    vtkDoubleArray* SMScalars = vtkDoubleArray::New();
    SMScalars->SetName("SpuriousFilter");
    SMScalars->SetNumberOfComponents(1);

	// Create a point set for the output
	vtkPoints * outputPoints = vtkPoints::New();
	output->SetPoints(outputPoints);
	outputPoints->Delete();

	// Create a line array for the output
	vtkCellArray * outputLines = vtkCellArray::New();
	output->SetLines(outputLines);
	outputLines->Delete();

	//
    //      PROCESS
    //

   // using namespace std;
    ofstream dataFile("/home/linux/Stephan/test.txt");

	// Number of points in the current fiber, and a list of its point IDs
	vtkIdType numberOfPoints;
	vtkIdType * pointList;

	// Setup progress bar
    int numberOfCells = 300;//inputLines->GetNumberOfCells();
    int progressStep = 5;//numberOfCells / 25;
    progressStep += (progressStep == 0) ? 1 : 0;

    double minDistSquared = ps->minDist * ps->minDist;

    double* fiberScores = ps->fiberScores;
    int* fiberNumPoints;
    int* fiberStartIds = ps->fiberStartIds;
    double* fiberData;
    double* fiberTangents;
    double* fiberMinScores = ps->fiberMinScores;
    double avgScoreTotal = ps->avgScoreTotal;



    // precompute the kernel
//    int numOrientations = 162;
//    int dims = 15;
//    double dimStep = 0.5;
//    int N = dims/dimStep;
//    this->SetProgressText("Precomputing the kernel...");
//    this->UpdateProgress(0.0);
//    int c = 0;
//    for(int i = 0; i<numOrientations; i++)
//    {
//        if ((i % 5) == 0)
//        {
//            this->UpdateProgress((double) i / (double) numOrientations);
//        }
//
//        double* arg2 = EulerAngles(&orientations[i]);
//
//        for(double x = 0; x<dims; x+=dimStep)
//        {
//            for(double y = 0; y<dims; y+=dimStep)
//            {
//                for(double z = 0; z<dims; z+=dimStep)
//                {
//                    double kernelval = kernel(x,
//                                              y,
//                                              z,
//                                              arg2[0],
//                                              arg2[1],
//                                              this->ps->D33, this->ps->D44, this->ps->t);
//                    kernelOutput[c] = kernelval;
//                    c++;
//
//                }
//            }
//        }
//        free(arg2);
//    }


    //double* x, double* y, double* r, double* v)



    this->SetProgressText("Filtering fibers...");
    this->UpdateProgress(0.0);

    size_t sizeOfDouble1 = sizeof(double);
    size_t sizeOfDouble3 = sizeof(double)*3;

	/** COMPUTE **/
	if(this->ps->requireRecompute)
	{
        // save fibers in raw form
        fiberNumPoints = (int*)malloc(numberOfCells*sizeof(int));
        if(fiberStartIds != NULL)
            free(fiberStartIds);
        fiberStartIds = (int*)malloc(numberOfCells*sizeof(int));
        int totalNumPoints = 0;
        for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
        {
            vtkCell * currentCell = input->GetCell(lineId);
            int numberOfFiberPoints = currentCell->GetNumberOfPoints();
            fiberNumPoints[lineId] = numberOfFiberPoints;

            fiberStartIds[lineId] = totalNumPoints*3;
            totalNumPoints+=numberOfFiberPoints;
        }

        // copy fiber data and compute tangents
        fiberData = (double*)malloc(totalNumPoints*sizeOfDouble3);   // x y z
        fiberTangents = (double*)malloc(totalNumPoints*sizeOfDouble3);   // dx dy dz
        for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
        {
            // Get the data of the current fiber
            vtkCell * currentCell = input->GetCell(lineId);

            // Get cell containing fiber
            inputLines->GetNextCell(numberOfPoints, pointList);
            int numberOfFiberPoints = currentCell->GetNumberOfPoints();

            // Previous point coordinates
            double* prev_p = (double*) malloc(sizeOfDouble3);

            // Current point coordinates
            double p[3];

            int fiberStartId = fiberStartIds[lineId];

            // Loop through all points in the fiber
            for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
            {
                // Get the point ID of the current fiber point
                vtkIdType currentPointId = currentCell->GetPointId(pointId);

                // Copy the point coordinates to the output
                inputPoints->GetPoint(currentPointId, p);

                // copy point to fast array
                fiberData[fiberStartId+pointId*3] = p[0];
                fiberData[fiberStartId+pointId*3+1] = p[1];
                fiberData[fiberStartId+pointId*3+2] = p[2];

                if(pointId > 0)
                {
                    double* tangent = Normalize(Difference(p,prev_p));
                    fiberTangents[fiberStartId+pointId*3] = tangent[0];
                    fiberTangents[fiberStartId+pointId*3+1] = tangent[1];
                    fiberTangents[fiberStartId+pointId*3+2] = tangent[2];
                }
                else
                {
                    fiberTangents[fiberStartId+pointId*3] = 0.0;
                    fiberTangents[fiberStartId+pointId*3+1] = 0.0;
                    fiberTangents[fiberStartId+pointId*3+2] = 0.0;
                }

                // Set previous point
                memcpy(prev_p,p,3*sizeof(double));
            }

            // release memory
            free(prev_p);
        }

        clock_t lastTime = clock();


        omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions

        // prepare cache to utilise kernel associative property ab=ba where a,b in R3/S2
        //int* cacheIndices = (int*)malloc(totalNumPoints*2*sizeof(int)); // number of unique connections
        double* cacheValues = (double*)malloc(totalNumPoints*totalNumPoints*sizeof(double)); // diagonally symmetric matrix of output values
        int sfac1 = 0,sfac2 = 0;
        printf("total num points: %d \n",totalNumPoints);

        // Compute kernel for fibers
        if(fiberScores != NULL)
            free(fiberScores);
        fiberScores = (double*)malloc(totalNumPoints*sizeOfDouble1);
        for (int lineId = 0; lineId < numberOfCells; ++lineId)
        {
            // Update the progress bar
            if ((lineId % progressStep) == 0)
            {
                this->UpdateProgress((double) lineId / (double) numberOfCells);
            }

            // Get the data of the current fiber
            int numberOfFiberPoints = fiberNumPoints[lineId];

            // Current point coordinates
            int fiberStartId = fiberStartIds[lineId];

            double* p;
            double* tangent;
            double* tangentFlipped;

            // Loop through all points in the fiber
            for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
            {
                double fscore = 0.0;
                if(pointId > 0)
                {
                    int id1 = fiberStartId + (pointId)*3;
                    int id11 = id1/3;

                    // load point
                    p = fiberData + id1;
                    // load f1 tangent
                    tangent = fiberTangents + id1;

                    // compute flipped f1 tangent
                    if(ps->applyKernelInBothDirs)
                        tangentFlipped = Flip(tangent);

                    // second loop around f2
                    int totalNumPoint = 0;
                    double* scoringResults = (double*)calloc(numberOfCells,sizeOfDouble1);
                    #pragma omp parallel for
                    for (int lineIdf2 = 0; lineIdf2 < numberOfCells; ++lineIdf2)
                    {
                        // skip equal fibers to prevent unnecessary calculations
                        if(lineId == lineIdf2)
                            continue;

                        double* pf2;
                        double* tangentf2;

                        int numberOfFiberPointsf2 = fiberNumPoints[lineIdf2];
                        int fiberStartIdf2 = fiberStartIds[lineIdf2];

                        double score = 0.0;
                        for (int pointIdf2 = 1; pointIdf2 < numberOfFiberPointsf2; ++pointIdf2)
                        {
                            int id2 = fiberStartIdf2 + pointIdf2*3;
                            int id22 = id2/3;
                            pf2 = fiberData + id2;
                            tangentf2 = fiberTangents + id2;

                            double sval = 0.0;
                            // is a cached result available?
//                            if(id22 < id11)
//                            {
//                                printf("request at: (%d,%d):%d of max:%dx%d=%d\n",id11,id22,id11*totalNumPoints+id22,totalNumPoints,totalNumPoints,totalNumPoints*totalNumPoints);
//                                sval = cacheValues[id11*totalNumPoints+id22];
//                                sfac1++;
//                            }
//
//                            // if not, compute the kernel
//                            else
//                            {
                                // check distance
                                if(Norm(Subtract3(p,pf2)) + 10.0*Norm(Subtract3(tangent, tangentf2)) < minDistSquared)
                                {
                                    // evaluate kernel
                                    //printf("forward: %f  backward:%f \n",k2(p,pf2,tangent,tangentf2),k2(pf2,p,tangent,tangentf2));
                                    sval = k2(p,pf2,tangent,tangentf2);

                                    // flipped kernel and take largest response
                                    if(ps->applyKernelInBothDirs)
                                    {
                                        double sval2 = k2(p,pf2,tangentFlipped,tangentf2);
                                        if(sval2 > sval)
                                            sval = sval2;
                                    }

                                    //printf("fill at: (%d,%d):%d of max:%dx%d=%d\n",id1,id2,id1*totalNumPoints+id2,totalNumPoints,totalNumPoints,totalNumPoints*totalNumPoints);
                                    cacheValues[id11*totalNumPoints+id22] = sval;
                                }
//                            }
//                            sfac2++;

                            // add scoring value. check against NaN
                            if(!(sval!=sval))
                                score += sval;
                        }
                        scoringResults[lineIdf2] = score;
                    }

                    // sum individual scorings
                    for(int i = 0; i<numberOfCells; i++)
                    {
                        //dataFile << scoringResults[i] << "\n";
                        fscore += scoringResults[i];
                    }
                    free(scoringResults);
                }

                fiberScores[fiberStartId/3 + pointId] = fscore;

            }

            // compute remaining time
            clock_t newTime = clock();
            QString progtext = QString("Filtering fibers... (time remaining = %1 seconds)").arg(int(double(newTime - lastTime) / CLOCKS_PER_SEC * (numberOfCells - lineId)));
            this->SetProgressText(progtext.toLocal8Bit().data());
            lastTime = newTime;

            if(this->GetAbortExecute())
                return;
        }

        free(cacheValues);
        printf("speedup factor: %f\n",(double)sfac1/sfac2);


        /** compute windows **/

        // window settings
        double meanNumPoints = 0.0;
        int minNumPoints = 1e9;
        for (int lineId = 0; lineId < numberOfCells; ++lineId)
        {
            int numberOfFiberPoints = fiberNumPoints[lineId];
            meanNumPoints += numberOfFiberPoints;
            if(numberOfFiberPoints < minNumPoints)
                minNumPoints = numberOfFiberPoints;
        }
        meanNumPoints /= numberOfCells;
        double windowSize = ps->windowSize;
        int adjustedWindowSize = std::min(minNumPoints,int(meanNumPoints/windowSize + .5)); // adjust window size

        // compute the windows
        avgScoreTotal = 0.0;
        if(fiberMinScores != NULL)
           free(fiberMinScores);
        fiberMinScores = (double*)malloc(numberOfCells*sizeof(double));
        for (int lineId = 0; lineId < numberOfCells; ++lineId)
        {
            int numberOfFiberPoints = fiberNumPoints[lineId];

            // Loop through all points in the fiber within the window
            // skip first point due to irrelevant score
            int windowStart = 1;

            // Current point coordinates
            int fiberStartId = fiberStartIds[lineId];

            double minScore = 1e99;
            //double avgScore = 0.0;

            while(true)
            {
                int windowEnd = windowStart + adjustedWindowSize - 1;
                if(windowEnd >= numberOfFiberPoints - 1)
                    break;

                double windowSum = 0.0;
                for (int pointId = windowStart; pointId <= windowEnd; ++pointId)
                {
                    double score = fiberScores[fiberStartId/3 + pointId];
                    windowSum += score;
                }

                windowSum /= windowEnd - windowStart + 1;

                //printf("windowsum: %f \n",windowSum);

                if(windowSum < minScore)
                    minScore = windowSum;
                //avgScore += windowSum;

                windowStart++;
            }

            //avgScore /= windowStart;
            //avgScoreTotal += avgScore;

            fiberMinScores[lineId] = minScore;
            printf("minscore: %f \n",minScore);
        }

        for (int pointId = 0; pointId < totalNumPoints; ++pointId)
        {
            avgScoreTotal += fiberScores[pointId];
            dataFile << fiberScores[pointId] << " " << "\n";
        }
        avgScoreTotal /= totalNumPoints;

        printf("intlength: %d, averageTotal: %f \n",adjustedWindowSize,avgScoreTotal);

        free(fiberData);
        free(fiberNumPoints);
        free(fiberTangents);

        ps->fiberScores = fiberScores;
        ps->fiberStartIds = fiberStartIds;
        ps->fiberMinScores = fiberMinScores;
        ps->avgScoreTotal = avgScoreTotal;

        // compute goodness
        ps->goodness = avgScoreTotal;
    }
    /** END COMPUTE **/

	// save fiber results
	double epsilon = ps->epsilon;
	printf("filter with epsilon: %f \n", epsilon);
	for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
	{
	    // Update the progress bar
		if ((lineId % progressStep) == 0)
		{
			this->UpdateProgress((double) lineId / (double) numberOfCells);
		}

		// threshold the score against epsilon
		//printf("sdaasdasd %f %f \n",fiberMinScores[lineId], ps->avgScoreTotal);
        if( fiberMinScores[lineId] / avgScoreTotal < epsilon)
            continue;
        //if( fiberMinScores[lineId] < epsilon)
        //    continue;

        // Get the data of the current fiber
        vtkCell * currentCell = input->GetCell(lineId);
        int numberOfFiberPoints = currentCell->GetNumberOfPoints();

		// Create an ID list for the output fiber
		vtkIdList * newFiberList = vtkIdList::New();

		// Current point coordinates
		double p[3];

        int fiberStartId = fiberStartIds[lineId];

        // Loop through all points in the fiber
		for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
		{
            // Get the point ID of the current fiber point
			vtkIdType currentPointId = currentCell->GetPointId(pointId);

			// Copy the point coordinates to the output
			inputPoints->GetPoint(currentPointId, p);

            vtkIdType newPointId = outputPoints->InsertNextPoint(p);
            newFiberList->InsertNextId(newPointId);

            // include old scalar values
            for(int i = 0; i < numberOfScalarTypes; i++)
            {
                // Get the scalar value
                double scalar = inputPD->GetArray(i)->GetTuple1(currentPointId);

                // Copy the scalar value to the output
                outputScalarsList.at(i)->InsertNextTuple1(scalar);
            }

            // include new scoring
            SMScalars->InsertNextTuple1(fiberScores[fiberStartId/3 + pointId]);
		}

		// Add the new fiber to the output
		outputLines->InsertNextCell(newFiberList);
	}

	// Add scalar arrays
	for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = outputScalarsList.at(i);
        output->GetPointData()->AddArray(outputScalars);
    }

    // Add SM scalar array
    output->GetPointData()->AddArray(SMScalars);
    output->GetPointData()->SetActiveScalars("SpuriousFilter");

	// Finalize the progress bar
	this->UpdateProgress(1.0);
    this->ps->requireRecompute = false;

    dataFile.close();
}

} // namespace bmia
