#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

//center of cordinate is on the centeral axis

const double number_of_primart_photon=1e6-1700; //mev
const double _total_mass_attenuation_coefficient_0_1MeV=0.167; //cm2/g
const double _total_mass_attenuation_coefficient_0_3MeV=0.118; //cm2/g
const double _total_mass_attenuation_coefficient_0_5MeV=0.0966; //cm2/g
const double _total_mass_attenuation_coefficient_1_25MeV=0.0630; //cm2/g
const double _total_mass_attenuation_coefficient_0_8MeV=0.0786; //cm2/g
const double _total_mass_attenuation_coefficient_1_0MeV=0.0706; //cm2/g
const double _total_mass_attenuation_coefficient_1_5MeV=0.0575; //cm2/g
const double _total_mass_attenuation_coefficient_2_0MeV=0.0493; //cm2/g
const double _total_mass_attenuation_coefficient_3_0MeV=0.0396; //cm2/g
const double _total_mass_attenuation_coefficient_4_0MeV=0.0339; //cm2/g
const double _total_mass_attenuation_coefficient_5_0MeV=0.0301; //cm2/g
const double fluence=20000; // photon/cm2
const double _phantom_dimention[]={20,20,20}; //cm
const bool is_circular_field_size=false;
const double _field_size_radius=5; //cm
const bool is_rectangular_field_size=true;
const double rectangular_field_size[]={10,10}; //cm
const double depth=15; //cm


double calculate_middle_of_the_voxel_r(double r_id)
{
        //calculate middle of the voxel, depth
    double R=0;
    if(r_id!=0)
    {
        double r0 = 0.05*(2*(r_id-1)*(r_id-1)+1); //cm
        double r1 = 0.05*(2*(r_id)*(r_id)+1); //cm
        R=(r0+r1)/2;
    }
    else
    {
        R=0.025; //cm
    }

    return R;
}

double calculate_middle_of_the_voxel_theta(double interaction_point_theta_id)
{
     //calculate middle of the voxel, theta
    double interaction_point_theta=0;
    if(interaction_point_theta_id!=0)
    {
        double theta0=(interaction_point_theta_id-1)*3.75;
        double tehta1=(interaction_point_theta_id)*3.75;
        interaction_point_theta=(theta0+tehta1)/2;
    }
    else
    {
            double interaction_point_theta=1.825;
    }

    return interaction_point_theta;   
}

TH2D* CalculateKernel()
{
        //read kernels
        TFile* f =new TFile("kernel_4_0mev.root");
        TH2D* edepPrimaryH = (TH2D*)f->Get("2d_PrimaryedepKernels");
        TH2D* edepFirstH = (TH2D*)f->Get("2D_FirstedepKernels");
        TH2D* edepSecondH = (TH2D*)f->Get("2D_SecondedepKernels");
        TH2D* edepMultipleH = (TH2D*)f->Get("2D_MultipleedepKernels");
        TH2D* edepeBremAnnihilH = (TH2D*)f->Get("2D_eBremAnnihiledepKernels");
        TH2D* massH = (TH2D*)f->Get("2d_Mass");

        TH2D* _result=new TH2D("kernels","kernels",25,0,25,48,0,48); 

        //mev/photon
        for(int j=0;j<25;j++)
            for(int i=0;i<48;i++)
                {
                    double edk=((edepPrimaryH->GetBinContent(i,j)+
                            edepFirstH->GetBinContent(i,j)+
                            edepSecondH->GetBinContent(i,j)+
                            edepMultipleH->GetBinContent(i,j)+
                            edepeBremAnnihilH->GetBinContent(i,j))/number_of_primart_photon);
                    _result->SetBinContent(j,i,edk); //(r,theta)
                }

        return _result;
}

double GetEDK(double r,double interaction_point_theta_id,TH2D* KERNEL)
{
    //calculate deposition_point_theta
    double interaction_point_theta=interaction_point_theta_id*3.75;
    double deposition_point_theta=180-interaction_point_theta;
    double deposition_point_theta_id=deposition_point_theta/3.75;

    //get total kernel
    TH2D* _totalKernel=KERNEL;
    
    //get edk
    double edk=_totalKernel->GetBinContent(r,deposition_point_theta_id);

    return edk;
}

double Calculate_Number_of_interaction(double r ,double interaction_point_theta_id,double mass_attenuation_coefficient)   // photons/g
{
    double R=calculate_middle_of_the_voxel_r(r);
    double interaction_point_theta=calculate_middle_of_the_voxel_theta(interaction_point_theta_id);

    double DEPTH= depth + R*TMath::Cos((TMath::Pi()*interaction_point_theta)/(180));
    double number_of_interaction=mass_attenuation_coefficient*fluence*TMath::Exp(-mass_attenuation_coefficient*DEPTH); // photons/g

    return number_of_interaction;
}

bool check_voxel_have_interacted(double r_id,double theta_i,double lateral_distance_dose_calculation_point)
{
    bool result=false;
    double theta=calculate_middle_of_the_voxel_theta(theta_i);
    double R=calculate_middle_of_the_voxel_r(r_id);

    //caculate off axis distance
    double field_boundary=rectangular_field_size[0]/2;
    double off_axis_distance_min   = -R*TMath::Sin((theta*TMath::Pi())/(180))+lateral_distance_dose_calculation_point;
    double off_axis_distance_plus  = R*TMath::Sin((theta*TMath::Pi())/(180))+lateral_distance_dose_calculation_point;  

    if(off_axis_distance_min>-field_boundary && off_axis_distance_plus<field_boundary)result=true;
    else result=false;

    return result;
}

TTree* Normalize(TTree* Tdose)
{
    //find max
    double max=0;
    double dose;
    Tdose->SetBranchAddress("edep",&dose);
    for(int i=0; i<Tdose->GetEntries(); i++)
    {
        Tdose->GetEntry(i);
        if(max<dose)max=dose;
    }

    //create new TTree
    double pdd;
    double z;
    TTree* Tpdd = new TTree("Ntuple1", "Edep");
    Tpdd->Branch("edep", &pdd, "d/D");
    Tpdd->Branch("X", &z, "z/D");

    //loot to dose TTree and fill pdd
    double d;
    double Z;
    Tdose->SetBranchAddress("edep",&d);
    Tdose->SetBranchAddress("X",&Z);
    for(int i=0; i<Tdose->GetEntries(); i++)
    {
        Tdose->GetEntry(i);
        pdd=(d/max)*100;
        z=Z;
        Tpdd->Fill();
    }

        return Tpdd;
}

bool check_if_in_the_phantom(double r_id,double theta_id,double lateral_distance_dose_calculation_point)
{
    bool result=false;

    //calculate middle of the voxel
    double theta=calculate_middle_of_the_voxel_theta(theta_id);
    double R=calculate_middle_of_the_voxel_r(r_id);

    //calculate off axis distance and zepth respect to point of interaction
    double off_axis_distance_plus_right=lateral_distance_dose_calculation_point+R*TMath::Sin((theta*TMath::Pi())/(180));
    double off_axis_distance_min_left=lateral_distance_dose_calculation_point-R*TMath::Sin((theta*TMath::Pi())/(180));
    double z=R*TMath::Cos((theta*TMath::Pi())/(180));

    if(lateral_distance_dose_calculation_point<0)
        if(  off_axis_distance_min_left>=((-1)*_phantom_dimention[0])/2 && (z>(-1)*depth && z<_phantom_dimention[2]-depth)  )result=true;
        else result=false;
    else
        if( off_axis_distance_plus_right<=(_phantom_dimention[0])/2  && (z>(-1)*depth && z<_phantom_dimention[2]-depth)  )result=true;
        else result=false;
    return result;
}

void Draw_PDD(double mass_attenuation_coefficient)
{
double point[80];
double step=_phantom_dimention[0]/80; //cm

// generate points
for(int i=1;i<=80;i++)point[i-1]=((-1)*_phantom_dimention[0]/2)+(i-1)*step;

TH2D* KERNEL=CalculateKernel();

double dose;
double z;
TTree* pdd = new TTree("Ntuple1", "Edep");
pdd->Branch("edep", &dose, "edep/D");
pdd->Branch("X", &z, "X/D");


for(double lateral_distance : point)
{
    for(int r_id=0;r_id<25;r_id++)
        for(int theta_id=0;theta_id<48;theta_id++)
        {
            if(check_if_in_the_phantom(r_id,theta_id,lateral_distance))
                if(check_voxel_have_interacted(r_id,theta_id,lateral_distance))
                {
                    double No_interaction=Calculate_Number_of_interaction(r_id,theta_id,mass_attenuation_coefficient);
                    double edk=GetEDK(r_id,theta_id,KERNEL);
                    // if(lateral_distance<-5 || lateral_distance>5)
                    //     dose+=No_interaction*edk*0.5;
                    // else  dose+=No_interaction*edk;
                    dose+=No_interaction*edk;

                }
        }
    z=lateral_distance;
    pdd->Fill();
    dose=0;
}


auto Tpdd=Normalize(pdd);

TCanvas* c1 = new TCanvas("c1", "  ");
TFile* f =new TFile("4mev.root");
TTree* tree = (TTree*)f->Get("Ntuple1");
auto com=Normalize(tree);
Tpdd->SetMarkerStyle(26);
com->SetMarkerStyle(25); 
Tpdd->Draw("edep:X","","P"); 
TH2F *htemp = (TH2F*)c1->GetPrimitive("htemp"); 
htemp->SetTitle("Lateral Profile ; off axis distance (cm) ; Relative Dose");
com->Draw("edep:(X/5)-10","","sameP");
gPad->Update();


auto legend = new TLegend(0.1,0.7,0.48,0.9);
legend->SetHeader("4Mev 10x10 cm^{2}_{}","C");
legend->AddEntry(Tpdd,"Convolution Calculation","P");
legend->AddEntry(com,"MC Geant4","p"); 
legend->SetBorderSize(0);
legend->Draw();

c1->cd();
c1->Update(); 
}

void LATERAL()
{
    Draw_PDD(_total_mass_attenuation_coefficient_4_0MeV);
}