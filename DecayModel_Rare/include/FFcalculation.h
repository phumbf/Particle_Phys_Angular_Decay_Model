#ifndef FFCALCULATION_H 
#define FFCALCULATION_H 1

#include "TVectorD.h"
#include "FormFactorBase.h" 

#include <array>

//Classes to determine the different form factors as a function of kinematic terms. The basis are taken from the 
//different papers seen in the form_factors folder
//
class BSZ15FFCalculation : virtual public VectorFormFactorBase { 

public: 
  BSZ15FFCalculation( );

  BSZ15FFCalculation( const BSZ15FFCalculation& other ) ; 
    
  BSZ15FFCalculation( const std::vector<double>& params );
  
  BSZ15FFCalculation( TVectorD& params );

  ~BSZ15FFCalculation( ); ///< Destructor

  VectorFormFactorBase* clone() const override { 
    return new BSZ15FFCalculation( *this );
  }
  
  double A1( const double qsq ) const override ;
  
  double A12( const double qsq ) const override ;
  
  double A0( const double qsq ) const override ;
    
  double A2( const double qsq ) const override ;
  
  double T1( const double qsq ) const override ;
  
  double T2( const double qsq ) const override ;
  
  double T3( const double qsq ) const override ;
  
  double T23( const double qsq ) const override ;
  
  double T3tilde( const double qsq ) const override ;
  
  double V( const double qsq ) const override ;

  void print() const override ;

  void setParams( const TVectorD& params ) override ;

  void setParams( const std::vector< double >& params ) override ;

  TVectorD getParams() const override ;

 protected:
  
  double zseries( const double qsq,
		  const double a0,
		  const double a1,
		  const double a2 ) const ;

  
  double pole( const double qsq, const double mass ) const;
    
 private:

  std::array< double, 21 > params_;   
};


class BZ04FFCalculation : virtual public VectorFormFactorBase { 

public: 
  /// Standard constructor
  BZ04FFCalculation( );

  BZ04FFCalculation( const BZ04FFCalculation& other ) ; 
    
  BZ04FFCalculation( const std::vector<double>& params );

  BZ04FFCalculation( TVectorD& params );

  ~BZ04FFCalculation( ); ///< Destructor

  VectorFormFactorBase* clone() const override { 
    return new BZ04FFCalculation( *this );
  }
    
  
  double A1( const double qsq ) const override ;
  
  double A12( const double qsq ) const override ;
  
  double A0( const double qsq ) const override;
    
  double A2( const double qsq ) const override;
  
  double T1( const double qsq ) const override;
  
  double T2( const double qsq ) const override;
  
  double T3( const double qsq ) const override;
  
  double T23( const double qsq ) const override;
  
  double T3tilde( const double qsq ) const override; 
  
  double V( const double qsq ) const override;

  void setParams( const TVectorD& params ) override;
  
  void setParams( const std::vector<double>& params ) override;

  void print() const override;

  TVectorD getParams() const override;

 protected:
  
  double pole( const double qsq, const double r, const double msq ) const; 

  double sqpole( const double qsq, const double r, const double msq ) const;
    
 private:
  
  std::array< double,  19 > params_ ;
};



class AAS07FFCalculation : public ScalarFormFactorBase { 

public: 
  /// Standard constructor
  AAS07FFCalculation( );

  AAS07FFCalculation( const AAS07FFCalculation& other ) ; 
    
  AAS07FFCalculation( const std::vector<double>& params );

  AAS07FFCalculation( TVectorD& params );

  ~AAS07FFCalculation( ); ///< Destructor

  ScalarFormFactorBase* clone() const override { 
    return new AAS07FFCalculation( *this );
  }
   
  double fplus( const double qsq ) const override;

  double fminus( const double qsq ) const override;

  double fzero( const double qsq ) const override;

  double fT( const double qsq ) const override;
  
  void setParams( const TVectorD& params ) override;
  
  void setParams( const std::vector<double>& params ) override;

  void print() const override;

  TVectorD getParams() const override;

 protected:
    
  std::array< double,  9 > params_ ;
};

class ESS18FFCalculation : public VectorFormFactorBase { 

public: 
  /// Standard constructor
  ESS18FFCalculation( );

  ESS18FFCalculation( const ESS18FFCalculation& other ) ; 
    
  ESS18FFCalculation( const std::vector<double>& params );

  ESS18FFCalculation( TVectorD& params );

  ~ESS18FFCalculation( ); ///< Destructor

  VectorFormFactorBase* clone() const override { 
    return new ESS18FFCalculation( *this );
  }

  double A1( const double qsq) const override;

  double A2( const double qsq) const override;
  
  double A12( const double qsq ) const override;

  double A0( const double qsq) const override;

  double V( const double qsq) const override;

  double T1( const double qsq) const override;

  double T2( const double qsq) const override;

  double T3( const double qsq) const  override;

  double T23( const double qsq) const override;

  double T3tilde( const double qsq ) const override; 
   
  void setParams( const TVectorD& params ) override;
  
  void setParams( const std::vector<double>& params ) override;

  void print() const override;

  TVectorD getParams() const override;

 protected:
    
  std::array< double,  21 > params_ ;
};

class YANG10FFCalculation : public VectorFormFactorBase { 

public: 
  /// Standard constructor
  YANG10FFCalculation( );

  YANG10FFCalculation( const YANG10FFCalculation& other ) ; 
    
  YANG10FFCalculation( const std::vector<double>& params );

  YANG10FFCalculation( TVectorD& params );

  ~YANG10FFCalculation( ); ///< Destructor

  VectorFormFactorBase* clone() const override { 
    return new YANG10FFCalculation( *this );
  }

  double A1( const double qsq) const override;

  double A2( const double qsq) const override;
  
  double A12( const double qsq ) const override;

  double A0( const double qsq) const override;

  double V( const double qsq) const override;

  double T1( const double qsq) const override;

  double T2( const double qsq) const override;

  double T3( const double qsq) const  override;

  double T23( const double qsq) const override;

  double T3tilde( const double qsq ) const override; 
   
  void setParams( const TVectorD& params ) override;
  
  void setParams( const std::vector<double>& params ) override;

  void print() const override;

  TVectorD getParams() const override;

 protected:
    
  std::array< double,  21 > params_ ;
};

class MS00FFCalculation : public VectorFormFactorBase { 

public: 
  /// Standard constructor
  MS00FFCalculation( );

  MS00FFCalculation( const MS00FFCalculation& other ) ; 
    
  MS00FFCalculation( const std::vector<double>& params );

  MS00FFCalculation( TVectorD& params );

  ~MS00FFCalculation( ); ///< Destructor

  VectorFormFactorBase* clone() const override { 
    return new MS00FFCalculation( *this );
  }

  double A1( const double qsq) const override;

  double A2( const double qsq) const override;
  
  double A12( const double qsq ) const override;

  double A0( const double qsq) const override;

  double V( const double qsq) const override;

  double T1( const double qsq) const override;

  double T2( const double qsq) const override;

  double T3( const double qsq) const  override;

  double T23( const double qsq) const override;

  double T3tilde( const double qsq ) const override; 
   
  void setParams( const TVectorD& params ) override;
  
  void setParams( const std::vector<double>& params ) override;

  void print() const override;

  TVectorD getParams() const override;

 protected:
    
  std::array< double,  21 > params_ ;
};

class ABHH99FFCalculation : public VectorFormFactorBase { 

public: 
  /// Standard constructor
  ABHH99FFCalculation( );

  ABHH99FFCalculation( const ABHH99FFCalculation& other ) ; 
    
  ABHH99FFCalculation( const std::vector<double>& params );

  ABHH99FFCalculation( TVectorD& params );

  ~ABHH99FFCalculation( ); ///< Destructor

  VectorFormFactorBase* clone() const override { 
    return new ABHH99FFCalculation( *this );
  }

  double A1( const double qsq) const override;

  double A2( const double qsq) const override;
  
  double A12( const double qsq ) const override;

  double A0( const double qsq) const override;

  double V( const double qsq) const override;

  double T1( const double qsq) const override;

  double T2( const double qsq) const override;

  double T3( const double qsq) const  override;

  double T23( const double qsq) const override;

  double T3tilde( const double qsq ) const override; 
   
  void setParams( const TVectorD& params ) override;
  
  void setParams( const std::vector<double>& params ) override;

  void print() const override;

  TVectorD getParams() const override;

 protected:
    
  std::array< double,  28 > params_ ;
};


class CFW10FFCalculation : public ScalarFormFactorBase { 

public: 
  /// Standard constructor
  CFW10FFCalculation( );

  CFW10FFCalculation( const CFW10FFCalculation& other ) ; 
    
  CFW10FFCalculation( const std::vector<double>& params );

  CFW10FFCalculation( TVectorD& params );

  ~CFW10FFCalculation( ); ///< Destructor

  ScalarFormFactorBase* clone() const override { 
    return new CFW10FFCalculation( *this );
  }
   
  double fplus( const double qsq ) const override;

  double fminus( const double qsq ) const override;

  double fzero( const double qsq ) const  override;

  double fT( const double qsq ) const override;
  
  void setParams( const TVectorD& params ) override;
  
  void setParams( const std::vector<double>& params ) override;

  void print() const override;

  TVectorD getParams() const override;

 protected:
    
  std::array< double,  9 > params_ ;
};

#endif // FFCALCULATION_H
