#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <list>
#include <deque>
#include <algorithm>
#include <math.h>
using namespace std;
class  isotope
{ //a single isotope mass and its abundance - need to be able to quicksort this
 public:isotope(double _mass, double _abundance);
   isotope()
   {;
   };
   double mass;
   double abundance;
};

typedef isotope peak;   //recycle the isotope definition to hold spectral peaks too
const double ElectronMass = 5.4857990943E-4;    //the electron rest mass, for correcting for charged molecules
class  BaseSpectrum     //an abstract base class corresponding to a spectrum, or peak list
{
 private: protected: public:virtual ~ BaseSpectrum();
   virtual int size() = 0;      //return number of peaks
   virtual void resize(int _npeaks) = 0;        //set number of peaks
   virtual double intensity(double mass) = 0;   //returns intensity as a function of mass
   virtual void push_back(peak & ob) = 0;       //add a peak to the spectrum
};

class  PeakList:public BaseSpectrum
                                //a simple peaklist object ie pairs of mass,intensity values
{ //essentially a wrapper for a deque
 private: protected:deque < peak > peaks;
 public:PeakList();
   PeakList(int _npeaks);
   int    size();
   void   resize(int _npeaks);
   virtual double intensity(double mass);
          peak & operator[] (int idx);  //returns a reference to the peak data
   //these wrap some of the other underlying deque functionality
   void   push_back(peak & ob); //add peak to the end
   void   push_front(peak & ob);        //add peak to start
   void   pop_back();   //delete peak from end
   void   pop_front();  //delete peak from front
};

class  GaussianList:public PeakList
                                //a spectrum consisting of peaks, with a specified gaussian width (fwhm)
{
 private: protected: public:double fwhm;
          GaussianList();
          GaussianList(int _npeaks);
   double intensity(double mass);
};

class  IntensityList:public BaseSpectrum
                                //a spectrum with intensity values uniformly sampled across a mass range
{
 private: protected:deque < double >peaks;
 public:double low, high;
   //intensity values will be uniformly distributed across mass range low to high
          IntensityList();
          IntensityList(int _npeaks);
   int    size();
   void   resize(int _npeaks);
   virtual double intensity(double mass);
   double &operator[] (int idx);
   void   push_back(peak & ob); //add peak to the end
   void   push_front(peak & ob);        //add peak to start
   void   pop_back();   //delete peak from end
   void   pop_front();  //delete peak from front
};

class  AtomIsoAbun      //the isotope abundances of an individual atom or group of atoms
{
 private: protected:int    niso;
   //must not allow this to be changed directly due to memory (de)allocation requirement
 public:       AtomIsoAbun();
          AtomIsoAbun(int _z, int _niso);
          AtomIsoAbun(const AtomIsoAbun & ob);
          virtual ~ AtomIsoAbun();
   int    z;    //atomic number
   char   symbol[4];    //atom symbol
   isotope *isotopes;   //the array of isotopes
   void   Setniso(int _niso);   //set the number of isotopes
   int    Getniso();
   void   Sort_abundance();     //sort the isotopes according to abundance
          AtomIsoAbun & operator=(AtomIsoAbun & ob);
   bool   operator<(const AtomIsoAbun & ob) const;
   bool   operator==(const AtomIsoAbun & ob) const;
   friend istream & operator>>(istream & stream, AtomIsoAbun & ob);
   friend ostream & operator<<(ostream & stream, AtomIsoAbun & ob);
};

class  AtomDistribution //describes a distribution of atoms
{
 public:double mass, abundance;
   //total mass, and fractional abundance    (nb mass is NOT automatically calculated - as this object does not know about the isotope masses)
          vector < int >number; //number of atoms of each isotope
   double CalculateMass(AtomIsoAbun * ob);      //calculate mass based on supplied isotope masses
   friend ostream & operator<<(ostream & stream, AtomDistribution & ob);
};

class  AtomEnsemble     //an ensemble of isotope distributions - for a single element (group) only
{
 private: protected: public:AtomEnsemble();
   AtomEnsemble(int _z, int _niso);
          AtomEnsemble(AtomIsoAbun & _atomdef);
   AtomIsoAbun atomdef; //the atom and its isotopes
          list < AtomDistribution > distribution;       //the different compositions
   friend ostream & operator<<(ostream & stream, AtomEnsemble & ob);
};

class  MolDistribution  //a molecule with a specified isotopic composition and charge for all element types
{ //not store the individual isotope masses here to avoid duplication
 private:      //object does NOT own the AtomDistribution objects
 protected: public:double mass, abundance;
   int    z;    //the charge on the molecule - a bit wasteful as all molecules in an ensemble are going to have the same charge (although that is not necessarily so)
          vector < AtomDistribution * >atomlist;        //vector of isotopic compositions for each element(group)
   double CalculateMass();      //calculate mass based on sums of masses stored for each element of atomlist
   bool   operator<(const MolDistribution & ob) const;
   friend ostream & operator<<(ostream & stream, MolDistribution & ob);
};

class  MolComposition
{ //composition of a molecule, in terms of atoms described by AtomIsoAbun objects
 private:      //ie CxHyOz
 protected: public:MolComposition();
   MolComposition(const MolComposition & ob);
          MolComposition(int _totalnat);
          virtual ~ MolComposition();
   int    totalnat;     //number of different atoms types(isotope labels are counted as different)
   AtomIsoAbun *atom;   //atom isotope abundances for each atom
   int   *nat;  //array containing number of atoms of each type
   int    z;    //charge
   void   Settotalnat(int _totalnat);   //set the total number of atoms
          MolComposition & operator=(MolComposition & ob);
   friend istream & operator>>(istream & stream, MolComposition & ob);
   friend ostream & operator<<(ostream & stream, MolComposition & ob);
};

class  MolEnsemble
{
 private: protected: public:virtual ~ MolEnsemble();
   MolEnsemble();
   MolEnsemble(MolComposition & comp, double thr);
          MolEnsemble(const MolEnsemble & ob);
   AtomEnsemble *elements;      //array (owned by this object) containing the atom ensembles - the atomlist pointers will point into at the distributions here
   int    totalnat;     //number of atom types - for naming consistency with MolComposition
   enum
   { mass, masstocharge } spectype;     //this controls whether the the calculated spectrum is of mass or mass-to-charge ratio
          list < MolDistribution > molecule;
   int    GenerateEnsemble(MolComposition & comp, double thr);  //generates an ensemble of molecules with different isotope distributions, with thr as the probability threshold
   void   CalculateMasses();    //calculate masses of all molecules on the list, from the atomdefs (also updates the total masses of each element)
   friend ostream & operator<<(ostream & stream, MolEnsemble & ob);
          MolEnsemble & operator=(MolEnsemble & ob);
   void   MakeSpectrum(BaseSpectrum * obPtr, double degen);     //make a mass spectrum in the supplied spectrum object, merge groups of peaks within degen units
};

class  IsoCalc  // collects together a table of elements/isotopes, a molecule composition, molecule ensemble and a peaklist
{
 private:double thr, degen;
   //full width at half max for gaussian peaks, and the accuracy threshold for distribution calculation, and parameter for peak merging
          deque < AtomIsoAbun > atomtable;      // a list containing all of the elements and their isotope distributions
   MolComposition molecule;     // the molecular composition
   MolEnsemble ensemble;        // ensemble of isotopically defined molecules
   GaussianList masspeaks;      // the peaklist
   bool   atomsloaded, molloaded;       // a pair of flags for the status of the atomtable (true=have a valid atom table) and molecular composition (true=ok composition supplied)
   bool   normalized;   // true if normalization of peaklist has been carried out
   bool   calculated;   // true if distribution has been calculated for current atoms and molecule
   bool   autocalc;     // true for autocalculation on changing parameters (ie if thr, z or composition is changed, this will automatically trigger a recalculation of the distribution)
 protected: public:       virtual ~ IsoCalc();
          IsoCalc();
          IsoCalc(char *filename);      //constructor that reads in an atom table from a file
   int    ReadAtomTable(const char *filename);  //read the atom table from a file, returns 0 on success
   int    SetComposition(char *formula);        //convert a string containing a composition to a MolComposition, returns 0 on success
   int    Calculate();  //calculate the isotope distribution, return 0 on success
   int    GetNPeaks(int &_npeaks);      //number of peaks in the peaklist, returns 0 on success
   int    Mass(int n, double &m);       //mass (m) of peak no. n, return 0 on success
   int    Abundance(int n, double &abun);       //abundance of peak n, return 0 on success
   int    Peak(int n, double &mass, double &abun);      //supplies both mass and abundance in one go (return 0 on success)
   int    Normalize();  // normalize the intensity so that most intense peak = 100%
   int    SetFwhm(const double _fwhm);  //setter and getter for the peak full width at half max, return 0 on success
   int    GetFwhm(double &_fwhm);
   int    SetDegen(const double _degen);
   int    GetDegen(double &_degen);
   int    SetThr(double _thr);  // setter and getter for isotope calculation threshold
   int    GetThr(double &_thr);
   int    SetCharge(const int _z);
   int    GetCharge(int &_z);
   int    SetMassToCharge(const bool _mtoz);    //if true, calculates mass-to-charge instead of mass
   int    SetAutoCalc(const bool _autocalc);    //setter and getter for the autocalc flag
   int    GetAutoCalc(bool & _autocalc);
   int    Intensity(double m, double &i);       //supplies intensity i at mass m calculated using gaussians
};
