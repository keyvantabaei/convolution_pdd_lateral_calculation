#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

//center of cordinate is at the point if caculated dose

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

TH2D* CalculateKernel(TString filename)
{
        //read kernels
        TFile* f =new TFile(filename);
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

double Calculate_Number_of_interaction(double depth /* cm */,double r ,double interaction_point_theta_id,double mass_attenuation_coefficient)   // photons/g
{
    double R=calculate_middle_of_the_voxel_r(r);
    double interaction_point_theta=calculate_middle_of_the_voxel_theta(interaction_point_theta_id);

    double DEPTH= depth + R*TMath::Cos((TMath::Pi()*interaction_point_theta)/(180));
    double number_of_interaction=mass_attenuation_coefficient*fluence*TMath::Exp(-mass_attenuation_coefficient*DEPTH); // photons/g

    return number_of_interaction;
}

bool check_voxel_have_interacted(double r_id,double theta_i)
{
    bool result=false;
    double theta=calculate_middle_of_the_voxel_theta(theta_i);
    double R=calculate_middle_of_the_voxel_r(r_id);
    // double z=R*TMath::Cos((theta*TMath::Pi())/(180));

    //caculate off axis distance
    double off_axis_distance=R*TMath::Sin((theta*TMath::Pi())/(180));
    if(is_circular_field_size)
        if(2*(off_axis_distance*off_axis_distance)<=(_field_size_radius*_field_size_radius))result=true;
        else result=false;
    else if(is_rectangular_field_size)
        if(off_axis_distance<(rectangular_field_size[0])/2 && off_axis_distance<(rectangular_field_size[0])/2)result=true;
        else result=false;
    else cout<<"---------->> field size have not chosen .";

    return result;
}

TTree* Normalize(TTree* Tdose)
{
    //find max
    double max=0;
    double dose;
    Tdose->SetBranchAddress("d",&dose);
    for(int i=0; i<Tdose->GetEntries(); i++)
    {
        Tdose->GetEntry(i);
        if(max<dose)max=dose;
    }

    //create new TTree
    double pdd;
    double z;
    TTree* Tpdd = new TTree("Ntuple1", "Edep");
    Tpdd->Branch("d", &pdd, "d/D");
    Tpdd->Branch("z", &z, "z/D");

    //loot to dose TTree and fill pdd
    double d;
    double Z;
    Tdose->SetBranchAddress("d",&d);
    Tdose->SetBranchAddress("z",&Z);
    for(int i=0; i<Tdose->GetEntries(); i++)
    {
        Tdose->GetEntry(i);
        pdd=(d/max)*100;
        z=Z;
        Tpdd->Fill();
    }

        return Tpdd;
}

bool check_if_in_the_phantom(double r_id,double theta_id,double depth )
{
    bool result=false;

    //calculate middle of the voxel
    double theta=calculate_middle_of_the_voxel_theta(theta_id);
    double R=calculate_middle_of_the_voxel_r(r_id);

    //calculate off axis distance and zepth respect to point of interaction
    double off_axis_distance=R*TMath::Sin((theta*TMath::Pi())/(180));
    double z=R*TMath::Cos((theta*TMath::Pi())/(180));

    if( off_axis_distance<=(_phantom_dimention[0])/2  && (z+depth>0 && z+depth<_phantom_dimention[2])  )result=true;
    else result=false;

    return result;
}

TTree* Get_TTree_PPD(TH2D* KERNEL,double mass_attenuation_coefficient)
{

// generate points
int point_number=100;
double point[point_number];
double step=_phantom_dimention[0]/point_number; //cm
for(int i=1;i<=point_number;i++)point[i-1]=(i-1)*step;

//Create TTree PDD
double dose;
double z;
TTree* pdd = new TTree("Ntuple1", "Edep");
pdd->Branch("d", &dose, "d/D");
pdd->Branch("z", &z, "z/D");

for(double depth : point)
{
    for(int r_id=0;r_id<25;r_id++)
        for(int theta_id=0;theta_id<48;theta_id++)
        {
            if(check_if_in_the_phantom(r_id,theta_id,depth))
                if(check_voxel_have_interacted(r_id,theta_id))
                {
                    double No_interaction=Calculate_Number_of_interaction(depth,r_id,theta_id,mass_attenuation_coefficient);
                    double edk=GetEDK(r_id,theta_id,KERNEL);
                    dose+=No_interaction*edk;

                }
        }
    z=depth;
    pdd->Fill();
    dose=0;
}

return pdd;
}

void Draw_PDD(double mass_attenuation_coefficient,double mass_attenuation_coefficient2,TString monte_carlo_data_name,TString monte_carlo_data_name2,TString kernel_data_name,TString kernel_data_name2)
{

//get kernels
TH2D* KERNEL=CalculateKernel(kernel_data_name);
TH2D* KERNEL2=CalculateKernel(kernel_data_name2);


auto Tpdd=Normalize(Get_TTree_PPD(KERNEL,mass_attenuation_coefficient));
auto Tpdd2=Normalize(Get_TTree_PPD(KERNEL2,mass_attenuation_coefficient2));

//create canvas
TCanvas* c1 = new TCanvas("c1", "  ");

//get monte carlo pdds
TFile* monte_carlo_data =new TFile(monte_carlo_data_name);
TTree* tree = (TTree*)monte_carlo_data->Get("Ntuple1");
TFile* monte_carlo_data2 =new TFile(monte_carlo_data_name2);
TTree* tree2 = (TTree*)monte_carlo_data2->Get("Ntuple1");

//monte carlo pdd to histograms
TH1D* hist=new TH1D("h","s",100,0,20);
double z;
double pdd;
tree->SetBranchAddress("edep",&pdd);
tree->SetBranchAddress("Z",&z);
for(int i=0; i<tree->GetEntries(); i++)
{
    tree->GetEntry(i);
    hist->SetBinContent(z,pdd);
}
//2
TH1D* hist2=new TH1D("h","s",100,0,20);
double z2;
double pdd2;
tree2->SetBranchAddress("edep",&pdd2);
tree2->SetBranchAddress("Z",&z2);
for(int i=0; i<tree2->GetEntries(); i++)
{
    tree2->GetEntry(i);
    hist2->SetBinContent(z2,pdd2);
}
//set dot style for data
Tpdd->SetMarkerStyle(21);
Tpdd2->SetMarkerStyle(22);
hist->SetMarkerStyle(25); 
hist2->SetMarkerStyle(26); 

//draw convolution calculation pdd
Tpdd->Draw("d:z","","P"); 
Tpdd2->Draw("d:z","","Psame"); 

//set title and axis lables
TH2F *htemp = (TH2F*)c1->GetPrimitive("htemp"); 
htemp->SetTitle("Percentage Depth Dose ; Depth (cm) ; Relative Dose");


//draw monte carlo pdd
hist->Draw("sameP");
hist2->Draw("sameP");
gPad->Update();


auto legend = new TLegend(0.1,0.7,0.48,0.9);
legend->SetHeader("","");
legend->AddEntry(Tpdd,"Convolution Calculation 4MeV","p");
legend->AddEntry(hist,"MC Geant4 4MeV","p"); 
legend->AddEntry(Tpdd2,"Convolution Calculation 1.25MeV","p");
legend->AddEntry(hist2,"MC Geant4 1.25MeV","p"); 
legend->SetBorderSize(0);
legend->Draw();

c1->cd();
c1->Update(); 
}

void PDD()
{
    //import kernel data
    TString kernel_data_name="kernel_4_0mev.root";
    TString kernel_data_name2="kernel_1_25mev.root";

    //import monte carlo data
    TString monte_carlo_data_name="pdd_4_0mev_10x10.root";
    TString monte_carlo_data_name2="pdd_1_25mev_10x10.root";

    Draw_PDD(_total_mass_attenuation_coefficient_4_0MeV,_total_mass_attenuation_coefficient_1_25MeV,monte_carlo_data_name,monte_carlo_data_name2,kernel_data_name,kernel_data_name2);
    // Draw_PDD(_total_mass_attenuation_coefficient_1_25MeV,monte_carlo_data_name2,kernel_data_name2);

}