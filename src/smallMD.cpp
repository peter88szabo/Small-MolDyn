// scope of local variables
// and global variables together
#include<iostream>
#include<cmath>
using std::cin;
using std::cout;
using std::endl;
using std::string;
 
// global variable
const int ndim=1;
int istep;
int maxstep=200;
int nfix=0;
double dt=0.1,dq=0.002;
double Etot,Ekin,Epot;
double Tini,Ttarg,tau;
double wmass[ndim];
double q[ndim];
double p[ndim];
double force[ndim];


void Velocity_Verlet();
void Force_calc();
void Vpotential();
void Kinetic_energy;
void Energy_calc();
void Init_parameters();
void Initialize_momenta(Tini);
void Thermostat_Berendsen();
void Act_temp(nfix);
  

int main()
{  
     Init_parameters();

     Energy_calc();
     cout << 0  << "    " <<q[0] << "  "<< p[0] << "  " << Epot << "   "<<Etot << endl;

     for(istep=1; istep=maxstep; istep++){
	Velocity_Verlet();
	Energy_calc();
        //cout << istep, q[0],p[0],Ekin,Epot,Etot << endl;
        cout << istep << "    " <<q[0] << "  "<< p[0] << "  " << Epot << "   "<<Etot << endl;
     }

}
///////////////////////////////////////////////////////////////////////////////

//==================================================================
void Init_parameters(){

    int i;
    for(i=0; i<ndim; i++){
       q[i]=2.0;
       p[i]=0.0;
       wmass[i]=1.0;
    }

}
//==================================================================

//==================================================================
void Force_calc(){
    int i;
    double Vp1,Vm1;

   // Numeric gradient 
    for(i=0; i<ndim; i++){

      //Vp1=V(q+dq)
        q[i] = q[i]+dq;
        Vpotential();

      //Vm1=V(q-dq)
        q[i] = q[i]-2.0*dq;
        Vpotential();

      //Force as the negative gradient
        force[i] = -0.5*(Vp1-Vm1)/dq;

      //restore q(i) coordinate (it has to be unchanged)
        q[i] = q[i]+dq;
    }

}
//==================================================================


//==================================================================
void Velocity_Verlet(){
     int i;

     Force_calc();

     for(i=0; i<ndim; i++){
        p[i] += 0.5*force[i]*dt;
        q[i] += p[i]/wmass[i]*dt;
     }

     Force_calc();

     for(i=0; i<ndim; i++){
        p[i] += 0.5*force[i]*dt;
     }

}
//==================================================================

//==================================================================
void Vpotential(){
     double rOH1,rOH2,theta;

     Cartesian_to_internal(rOH1,rOH2,theta);
     Epot = water_FF(rOH1,rOH2,theta);

}
//==================================================================

//==================================================================
void Kinetic_energy(){
    int i;
    double dum;

    dum=0;
    for(i=0; i<ndim; i++){
       dum += p[i]*p[i]/wmass[i]/2.0;
    }

    Ekin=dum;

}
//==================================================================



//==================================================================
void Energy_calc(){

    Vpotential();
    Kinetic_energy();

    Etot=Ekin+Epot;
}
//==================================================================
//


//=================================================================
//     Berendesen thermostat:
//     rescaling of momenta to get the
//     correct momenta at the desired temperature
//=================================================================
void Thermostat_Berendsen(){
   
    int i;	  
    double Tact,scal;

  //calculate the acutal temperature:
    Tact=Act_temp();

  //scaling factor
    scal = sqrt(1.0 + dt*(Ttarg-Tact)/Tact/tau);

 // scaling of momenta (velocities)
    for(i=0; i<ndim; i++){
       p[i] = p[i]*scal;
    }

}
//=================================================================


//==================================================================
//      Initialize atomic momenta to fulfill
//      the Maxwell-Boltzmann distribution
//
//      The formula for velocities:
//      v(i)=sqrt(RT/m(i))*N(0,1)
//
//      where N(0,1) is a random number with normal distirbution
//
//      But we calculate momenta instead of velocities
//==================================================================
void Initialize_momenta(Tini){
    int i;
    double Tini,RT,rnd;

    RT=Rgas*Tini

    for(i=0; i<ndim; i++){

       rnd=random_gaussian()

       p[i]=sqrt(wmass[i]*RT)*rnd
    }
}
//==================================================================


//==================================================================
//     Actual temperature of the system
//
//     N*R*T/2 = Ekin
//
//     where N is the number of degrees of freedom of the system
//==================================================================
//     nfix is used for the number of the fixeddegrees of freedom
//     (e.g. nfix=3, if translational modes are frozen)
//==================================================================
void Act_temp(nfix){
    int nfix;

     Kinetic_energy();

     Tact=2.d0*Ekin/float(ndim-nfix)/Rgas
     
}
//=================================================================


//==================================================================
//   Conversion from Cartesian(q) to Internal(rOH1,rOH2,theta)
//   coordinates
//==================================================================
double Cartesian_to_internal()
     double x1,y1,z1;
     double x2,y2,z2;
     double rOH1,rOH2,theta,dot_hoh;

//----------------------------------------
//   This convention should be followed
//----------------------------------------
//   Order of atoms:
//   Atom 1. = O  [q(0-2)]
//   Atom 2. = H1 [q(3-5)]
//   Atom 3. = H2 [q(6-8)]
//----------------------------------------

  //[O-H1] relative distance vector:
    x1 = q[0]-q[3];
    y1 = q[1]-q[4];
    z1 = q[2]-q[5];

  //[O-H2] relative distance vector:
    x2 = q(0)-q(6);
    y2 = q(1)-q(7);
    z2 = q(2)-q(8);


   //O-H1 and O-H2 atomic distances
    rOH1 = sqrt(x1*x1 + y1*y1 + z1*z1)
    rOH2 = sqrt(x2*x2 + y2*y2 + z2*z2)

   //dot product of [O-H1] and [O-H2] relative distance vectors
     dot_hoh = x1*x2 + y1*y2 + z1*z2

   //angle between [O-H1] and [O-H2]
     theta=acos(dot_hoh/rOH1/rOH2) //in radian

   return (rOH1,rOH2,theta);
}
//==================================================================


//==================================================================
//     A simple force field for the water molecule
//------------------------------------------------------------------
//     All variables (rOH2, rOH2, Vpes) are in
//     atomic units (bohr, Hartree), theta in radian
//------------------------------------------------------------------
//     Ref:
//     J. Chem. Phys. 135, 224516 (2011)
//     doi: 10.1063/1.3663219
//==================================================================
double water_FF(rOH1,rOH2,theta){

    const double pi=3.1415926535897932384626433832795;	
    double rOH1,rOH2,theta,Vpes;
    double x,VOH1,VOH2,VHOH;

   //-------------------------------------------------------------
   // Fitting parameters of the force field
   //-------------------------------------------------------------
     const double Dr=432.581; //kJ/mol
     const double req=0.9419/0.5291772;  // Å-->bohr
     const double beta=2.287*0.5291772; // 1/Å --> 1/bohr
     const double Teq=107.4*pi/180.0;   // deg --> rad
     const double Ct=367.81; // kJ/mol/rad^2
   //-------------------------------------------------------------

   // Morse for O-H1 stretching
     x = 1.0-exp(-beta*(rOH1-req));
     VOH1 = Dr*x*x;

   // Morse for O-H2 stretching
     x = 1.d0-exp(-beta*(rOH2-req));
     VOH2 = Dr*x*x;

   // Harmonic bending for H-O-H angle
     VHOH = 0.5*Ct*(theta-Teq)*(theta-Teq);

   // Total potential energy
     Vpes = VOH1+VOH2+VHOH; //in kJ/mol
     Vpes = Vpes/2625.5; // Hartree

     return Vpes;
}


