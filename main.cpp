/*------------------------------
2018/01/03/
DE-rand-to-best-1-bin
動作させるためにはEigen3.3.4が必要 (http://eigen.tuxfamily.org/index.php?title=News:Eigen_3.3.4_released!)

Eigenディレクトリを以下のディレクトリに設置
rand-to-best-1-bin/
                ├ Eigen/
                └ main.cpp <- this file

参考文献: パラメータの相関を考慮した適応型差分進化アルゴリズムJADEの改良
( http://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/1939-14.pdf )

実行前に設定する項目
・各種パラメータ( population, scaling factor, CR , dimension , vector_size, generation_limit)
  - population          : 個体数
  - Scalling_factor     : F,F2とあるがFは大域的探索の強度，F2は局所的探索強度を 0~1の範囲で決定する
  - CR                  : 交叉率
  - dimension           : 次元数(vectorの個数)
  - vector_size         : vectorの要素数
  - generation_limit    : 世代数制限(最大繰り返し回数)

・class Differential_Evolution内のpublicな関数
  - void init()               : 制約条件を適用した初期値の決定
  - void evaluate()           : 最大化の場合は符号を逆に書き換える
  - void evaluate_of_best_x() : 最大化の場合は符号を逆に書き換える
  - double objective_func()   : 設定したい目的関数を決定( 評価値(double型)をreturn )

その他
・一個体は dimension x vector_sizeの行列になる
・他のDE戦略にしたい場合，DE_main関数内の以下の式を変更するとok
    x_c[k]=x[a]+ F2*(x_best-x[a])+ F*(x[b]-x[c]);
・終了条件は繰り返し回数によって決定される
・defaultのobjective_funcはSphere関数になっている
*///---------------------------


#include<iostream>
#include<vector>
#include<random>
#include"Eigen/Dense"


//-----------------各種パラメータ----------------------
const int population =100; //個体数
const double F = 0.8; //scaling factor, 差分ベクトル
const double F2 = 0.2; //scaling factor, 差分ベクトル
const double CR =0.95; // 交叉率
const int dimension = 16; //次元数(vectorの個数)
const int vector_size = 4; //vectorの要素数
const int generation_limit = 1000; //世代数制限


//----------------ここまで-----------------------------

//個体xのmatrix(16x4)の予約
typedef Eigen::Matrix<double, dimension,vector_size> individual;


class Differential_Evolution{
private:
  individual x[population];
  individual x_c[population];
  individual x_best;

  std::mt19937 mt{std::random_device{}()};

public:
  double objective_func(const individual & m){
    //Sphere_func
    double ans=1;
    int row=m.rows(),col=m.cols();
    for(int i=0;i<row;i++){
      for(int j=0;j<col;j++){
        ans+=m(i,j)*m(i,j);
      }
    }
    return ans;
  }
  void init(){
    std::uniform_real_distribution<double> u(-5,5);
    for(int i=0;i<population;i++){
      for(int j=0;j<dimension;j++){
        for(int k=0;k<vector_size;k++) x[i](j,k)= u(mt);
      }
    }
  }

  void evaluate_of_best_x(){
    double min=objective_func(x[0]);
    for(int i=0;i<population;i++){
      if(objective_func(x[i])<=min){
        min=objective_func(x[i]);
        x_best=x[i];
      }
    }
  }
  bool evaluate(int i){
    //return objective_func(x_c[i]) >= objective_func(x[i]); //最大化
    return objective_func(x_c[i]) <= objective_func(x[i]); //最小化
  }



  void DE_main(){
    std::uniform_int_distribution<int> choice(0,population-1);
    std::uniform_real_distribution<double> u(0.0,1.0);
    for(int i=0;i<population;i++){
      int a=-1,b=-1,c=-1;
      while(a==b || a==c || b==c){
        a=choice(mt); b=choice(mt); c=choice(mt);
      }
      int j=choice(mt);
      for(int k=0;k<population;k++){
        if(k==j || u(mt) < CR){
          x_c[k]=x[a]+ F2*(x_best-x[a])+ F*(x[b]-x[c]);// rand_to_best_1_bin
        }else{
          x_c[k]=x[k];
        }
      }
      if(evaluate(i)) x[i]=x_c[i];
    }
  }

  double best_ans(){  return objective_func(x_best); }

  void cout_best(){   std::cout<<x_best<<std::endl;   }


};

int main(void){
  Differential_Evolution DE;
  DE.init();
  DE.evaluate_of_best_x();

  for(int i=0;i<generation_limit;i++){
    std::cout<<i<<" ";
    DE.DE_main();
    DE.evaluate_of_best_x();
    std::cout<<DE.best_ans()<<std::endl;
  }

  std::cout<<generation_limit<<" "<<DE.best_ans()<<std::endl;
  DE.cout_best();

}
