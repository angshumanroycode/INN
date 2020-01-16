#include <Rcpp.h>
using namespace Rcpp;

double Norm(NumericVector x){
return(sqrt(sum(x*x)));
}

double fC(double x,int n){
  double y=(2*x-n-1)/n;
  int ySign;
  if(y>=0){
    ySign=1;
  }else{
    ySign=-1;
  }
  return -ySign*log(1-y*ySign);
}

// [[Rcpp::export]]
NumericVector INNC(List dlist,int n,int d){
  Function orderC("order");
  List L=List::create();
  for(int u=0;u<n;u++){
    NumericMatrix R(n-1,d);
    NumericMatrix O(n-1,d);
    for(int v=0;v<d;v++){
      NumericMatrix dlistv=dlist[v];
      NumericVector vec=dlistv(_,u);
      vec.erase(u);
      NumericVector vec2=orderC(vec);
      O(_,v)=vec2-1;
      NumericVector vec3=orderC(vec2);
      R(_,v)=vec3;
    }
    NumericMatrix M(d,d-1);
    for(int v=0;v<d;v++){
      NumericMatrix U(n-1,d-1);
      int wv=0;
      NumericVector vec=O(_,v);
      for(int w=0;w<d;w++)
        if(w!=v){
          NumericVector vec2=R(_,w);
          NumericVector vec3=vec2[vec];
          U(_,wv)=vec3;
          wv++;
        }
      NumericMatrix V=U(Range(0,n-3),_);
      for(int t=1;t<n-2;t++)
        for(int a=t;a<n-2;a++)
          for(int b=0;b<d-1;b++)
            if(U(a,b)>U(t-1,b))
              V(a,b)=V(a,b)-1;
      NumericVector uvsum(d-1);
      for(int a=0;a<n-2;a++)
        for(int b=0;b<d-1;b++)
          uvsum[b]=uvsum[b]+fC(V(a,b),(n-a-1));
      M(v,_)=uvsum;
    }
    L.push_back(M);
  }
  NumericMatrix LMat(d,d-1);
  for(int u=0;u<n;u++){
    NumericMatrix Mat=L[u];
    for(int v=0;v<d;v++)
      for(int w=0;w<d-1;w++)
        LMat(v,w)=LMat(v,w)+Mat(v,w);
  }
  NumericVector T(d);
  for(int v=0;v<d;v++){
    NumericVector Lvec=LMat(v,_);
    double nsmval=Norm(Lvec);
    T[v]=nsmval;
  }
  double sm=sum(T);
  double mx=max(T);
  NumericVector Tn=NumericVector::create(sm,mx);
  return Tn;
}

NumericVector SignM(NumericVector z){
  double nz=Norm(z);
  if(nz==0){
    int l=z.size();
    NumericVector zero(l);
    return zero;
  }else{
    return z/nz;
  }
}

NumericVector fS(NumericVector x){
  double nx=Norm(x);
  double lx=log(1-nx);
  NumericVector s=SignM(x);
  return s*lx;
}

// [[Rcpp::export]]
NumericVector INNS(List dlist,int n,int d){
  Function orderC("order");
  NumericMatrix Y(n-1,d);
  NumericMatrix O(n-1,d-1);
  List L=List::create();
  for(int u=0;u<n;u++){
    for(int v=0;v<d;v++){
      NumericMatrix Mat=dlist[v];
      NumericVector vec=Mat(_,u);
      vec.erase(u);
      Y(_,v)=vec;
    }
    NumericMatrix M(d,d-1);
    for(int v=0;v<d;v++){
      NumericVector vec=Y(_,v);
      NumericVector vec2=orderC(vec);
      for(int w=0;w<d;w++){
        if(w<v){
          NumericVector vec3=Y(_,w);
          NumericVector vec4=vec3[vec2-1];
          O(_,w)=vec4;
        }else if(w>v){
          NumericVector vec3=Y(_,w);
          NumericVector vec4=vec3[vec2-1];
          O(_,w-1)=vec4;
        }
      }
      NumericVector val(d-1);
      for(int s=0;s<n-2;s++){
        NumericVector vec3=O(s,_);
        NumericVector r(d-1);
        for(int t=s+1;t<n-1;t++){
          NumericVector vec4=O(t,_);
          NumericVector s=SignM(vec3-vec4);
          r=r+s;
        }
        r=r/(n-s-1);
        NumericVector tr=fS(r);
        val=val+tr;
      }
      M(v,_)=val;
    }
    L.push_back(M);
  }
  NumericMatrix LMat(d,d-1);
  for(int u=0;u<n;u++){
    NumericMatrix Mat=L[u];
    for(int v=0;v<d;v++)
      for(int w=0;w<d-1;w++)
        LMat(v,w)=LMat(v,w)+Mat(v,w);
  }
  NumericVector T(d);
  for(int v=0;v<d;v++){
    NumericVector Lvec=LMat(v,_);
    double nsmval=Norm(Lvec);
    T[v]=nsmval;
  }
  double sm=sum(T);
  double mx=max(T);
  NumericVector Tn=NumericVector::create(sm,mx);
  return Tn;
}

// [[Rcpp::export]]
List INNCtest(List dlist,int n,int d,int REP){
  Function sampleC("sample.int");
  NumericMatrix permval(REP,2);
  for(int rep=0;rep<REP;rep++){
    List dlistR=List::create();
    for(int i=0;i<d;i++){
      NumericMatrix D=dlist[i];
      IntegerVector ind=sampleC(n);
      ind=ind-1;
      NumericMatrix E(n,n);
      for(int u=0;u<n-1;u++)
        for(int v=u+1;v<n;v++){
          E(u,v)=D(ind(u),ind(v));
          E(v,u)=E(u,v);
        }
      dlistR.push_back(E);
    }
    permval(rep,_)=INNC(dlistR,n,d);
  }
  NumericVector val=INNC(dlist,n,d);
  NumericVector pval(2);
  for(int rep=0;rep<REP;rep++)
    for(int i=0;i<2;i++)
    if(val[i]<permval(rep,i))
      pval[i]=pval[i]+1;
  pval=(pval+1)/(REP+1.0);
  List L=List::create();
  L.push_back(permval);
  L.push_back(val);
  L.push_back(pval);
  return L;
}

// [[Rcpp::export]]
List INNStest(List dlist,int n,int d,int REP){
  Function sampleC("sample.int");
  NumericMatrix permval(REP,2);
  for(int rep=0;rep<REP;rep++){
    List dlistR=List::create();
    for(int i=0;i<d;i++){
      NumericMatrix D=dlist[i];
      IntegerVector ind=sampleC(n);
      ind=ind-1;
      NumericMatrix E(n,n);
      for(int u=0;u<n-1;u++)
        for(int v=u+1;v<n;v++){
          E(u,v)=D(ind(u),ind(v));
          E(v,u)=E(u,v);
        }
      dlistR.push_back(E);
    }
    permval(rep,_)=INNS(dlistR,n,d);
  }
  NumericVector val=INNS(dlist,n,d);
  NumericVector pval(2);
  for(int rep=0;rep<REP;rep++)
    for(int i=0;i<2;i++)
      if(val[i]<permval(rep,i))
        pval[i]=pval[i]+1;
  pval=(pval+1)/(REP+1.0);
  List L=List::create();
  L.push_back(permval);
  L.push_back(val);
  L.push_back(pval);
  return L;
}

/*
// [[Rcpp::export]]
NumericVector INN2(List dlist1,List dlist2,int n,int d){
  Function orderC("order");
  NumericVector T2stat(2*d);
  for(int v=0;v<d;v++){
    NumericMatrix dlistv1=dlist1[v];
    NumericMatrix dlistv2=dlist2[v];
    NumericVector vsum(2);
    for(int u=0;u<n;u++){
      NumericVector vec11=dlistv1(_,u);
      vec11.erase(u);
      NumericVector vec21=orderC(vec11);
      NumericVector o1=vec21-1;
      NumericVector r1=orderC(vec21);
      //
      NumericVector vec12=dlistv2(_,u);
      vec12.erase(u);
      NumericVector vec22=orderC(vec12);
      NumericVector o2=vec22-1;
      NumericVector r2=orderC(vec22);
      //
      NumericMatrix U(n-1,2);
      NumericVector r2o1=r2[o1];
      U(_,0)=r2o1;
      NumericVector r1o2=r1[o2];
      U(_,1)=r1o2;
      NumericMatrix V=U(Range(0,n-3),_);
      for(int t=1;t<n-2;t++)
        for(int a=t;a<n-2;a++)
          for(int b=0;b<2;b++)
            if(U(a,b)>U(t-1,b))
              V(a,b)=V(a,b)-1;
      for(int a=0;a<n-2;a++)
        for(int b=0;b<2;b++)
          vsum[b]=vsum[b]+fC(V(a,b),(n-a-1));
    }
    for(int b=0;b<2;b++){
      if(vsum[b]<0)
        vsum[b]=-vsum[b];
      T2stat[2*v+b]=vsum[b];
    }
  }
  double sm=sum(T2stat);
  sm=sm/2;
  double mx=max(T2stat);
  NumericVector Tn=NumericVector::create(sm,mx);
  return Tn;
}
*/

// [[Rcpp::export]]
NumericVector INN2(List dlist,List dlistm,int n,int d){
  Function orderC("order");
  NumericVector T2stat(d);
  for(int v=0;v<d;v++){
    NumericMatrix DV=dlist[v];
    NumericMatrix DMV=dlistm[v];
    double vsum=0;
    for(int u=0;u<n;u++){
      NumericVector vec1=DV(_,u);
      vec1.erase(u);
      NumericVector vec2=orderC(vec1);
      NumericVector o=vec2-1;
      //
      NumericVector vecm1=DMV(_,u);
      vecm1.erase(u);
      NumericVector vecm2=orderC(vecm1);
      NumericVector r=orderC(vecm2);
      //
      NumericVector U=r[o];
      NumericVector V=U[Range(0,n-3)];
      for(int t=1;t<n-2;t++)
        for(int a=t;a<n-2;a++)
          if(U[a]>U[t-1])
            V[a]=V[a]-1;
      for(int a=0;a<n-2;a++)
        vsum=vsum+fC(V[a],(n-a-1));
    }
    if(vsum<0)
      vsum=-vsum;
    T2stat[v]=vsum;
  }
  double sm=sum(T2stat);
  double mx=max(T2stat);
  NumericVector Tn=NumericVector::create(sm,mx);
  return Tn;
}

// [[Rcpp::export]]
List INN2test(List dlist,int n,int d,int REP){
  Function sampleC("sample.int");
  Function inn2C("inn2c");
  NumericMatrix permval(REP,2);
  for(int rep=0;rep<REP;rep++){
    List dlistR=List::create();
    for(int i=0;i<d;i++){
      NumericMatrix D=dlist[i];
      IntegerVector ind=sampleC(n);
      ind=ind-1;
      NumericMatrix E(n,n);
      for(int u=0;u<n-1;u++)
        for(int v=u+1;v<n;v++){
          E(u,v)=D(ind(u),ind(v));
          E(v,u)=E(u,v);
        }
      dlistR.push_back(E);
    }
    NumericVector vec=inn2C(dlistR,n,d);
    permval(rep,_)=vec;
  }
  NumericVector val=inn2C(dlist,n,d);
  NumericVector pval(2);
  for(int rep=0;rep<REP;rep++)
    for(int i=0;i<2;i++)
      if(val[i]<permval(rep,i))
        pval[i]=pval[i]+1;
  pval=(pval+1)/(REP+1.0);
  List L=List::create();
  L.push_back(permval);
  L.push_back(val);
  L.push_back(pval);
  return L;
}