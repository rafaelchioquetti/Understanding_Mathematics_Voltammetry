#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

double convolution(const vector<double>& current, int i, double h)
{
    double sum = 0.0;

    for(int j = 1; j < i; j++)
    {
        double term1 = (2.0/3.0)*(current[j-1]-current[j])*
        (pow(i-j+1,1.5)-pow(i-j,1.5));

        double term2 = 2*(current[j]*(i-j+1)-current[j-1]*(i-j))*
        (sqrt(i-j+1)-sqrt(i-j));

        sum += term1 + term2;
    }

    return sqrt(h/M_PI)*sum;
}

int main()
{

    const double h = 0.02;
    const int n = 5000;
    const double a = (4.0/3.0)*sqrt(h/M_PI);

    vector<double> current(2*n,0.0);
    vector<double> potential(2*n,0.0);
    vector<double> conv_prev(2*n,0.0);

    double g = 50;
    double g_min = -50;



    bool forward_scan = true;

    int choice;

    cout<<"Choose simulation (1-6): ";
    cin>>choice;

    if(choice == 4)   // PTET non-buffered example
    {

        double L, alpha;

        cout<<"Lambda: ";
        cin>>L;

        cout<<"alpha: ";
        cin>>alpha;

        double psi = 0;
        double conv = 0;

        for(int i=0;i<2*n;i++)
        {

            if(i==0)
            {
                psi = 0;
                conv = 0;
                potential[i] = g;
                
            }
            else
            {

                if(forward_scan)
                {
                    g -= h;

                    if(g <= g_min)
                        forward_scan = false;
                }
                else
                {
                    g += h;
                }

                conv = convolution(current,i,h);

                double f = (a/2.0)*psi + conv;

                psi = -(L * (f * exp(g) + f - 1)) /
                (a * exp(g) * L + a * L + exp(alpha * g));

               current[i] = psi;
                potential[i] = g;
                conv_prev[i] = f;
            }
        }
    }

        if(choice == 5)   // PTET non-buffered example
    {

        double L, alpha;

        cout<<"Lambda: ";
        cin>>L;

        cout<<"alpha: ";
        cin>>alpha;

        double psi = 0;
        double conv = 0;

        for(int i=0;i<2*n;i++)
        {

            if(i==0)
            {
                psi = 0;
                conv = 0;
                potential[i] = g;
                
            }
            else
            {

                if(forward_scan)
                {
                    g -= h;

                    if(g <= g_min)
                        forward_scan = false;
                }
                else
                {
                    g += h;
                }

                conv = convolution(current,i,h);

                double f = (a/2.0)*psi + conv;

                psi = (1-f) / (exp(g) + a);

                current[i] = psi;
                potential[i] = g;
                conv_prev[i] = f;
            }
        }
    }

        if(choice == 6)   // PTET non-buffered example
    {

        double L, alpha;

        cout<<"Lambda: ";
        cin>>L;

        cout<<"alpha: ";
        cin>>alpha;

        double psi = 0;
        double conv = 0;

        for(int i=0;i<2*n;i++)
        {

            if(i==0)
            {
                psi = 0;
                conv = 0;
                potential[i] = g;
                
            }
            else
            {

                if(forward_scan)
                {
                    g -= h;

                    if(g <= g_min)
                        forward_scan = false;
                }
                else
                {
                    g += h;
                }

                conv = convolution(current,i,h);

                double f = (a/2.0)*psi + conv;

                psi = (L - f * L) / (a * L + alpha * exp(g));;

                current[i] = psi;
                potential[i] = g;
                conv_prev[i] = f;
            }
        }
    }

            if(choice == 7)   // PTET non-buffered example
    {

        double L, alpha;

        cout<<"Lambda: ";
        cin>>L;

        cout<<"alpha: ";
        cin>>alpha;

        double psi = 0;
        double conv = 0;

        for(int i=0;i<2*n;i++)
        {

            if(i==0)
            {
                psi = 0;
                conv = 0;
                potential[i] = g;
                
            }
            else
            {

                if(forward_scan)
                {
                    g -= h;

                    if(g <= g_min)
                        forward_scan = false;
                }
                else
                {
                    g += h;
                }

                conv = convolution(current,i,h);

                double f = (a/2.0)*psi + conv;

                psi = (L - f * L) / (a * L + alpha * exp(alpha * g));

                current[i] = psi;
                potential[i] = g;
                conv_prev[i] = f;
            }
        }
    }


        ofstream file("voltammogram.csv");

        file<<"Potential,Current\n";

        for(int i=0;i<2*n;i++)
            file<<potential[i]<<"\t"<<current[i]<<"\n";

        file.close();

        cout<<"Simulation finished. Data saved to voltammogram.csv\n";
    

}