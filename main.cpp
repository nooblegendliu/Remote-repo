#include <iostream>
#include <cmath>

#define u_1 0.0
#define u_2 0.0
#define rho_1 1.0
#define rho_2 0.125
#define P_1 1.0
#define P_2 0.1
#define gamma 1.4


class Solve_System{
private:
    // needed parameters
    const double c_1 = sqrt(gamma*P_1/rho_1);
    const double c_2 = sqrt(gamma*P_2/rho_2);
    const double alpha_1_init = rho_1 * c_1;
    const double alpha_2_init = rho_2 * c_2;
    // contact interruption
    double Pstar;
    double Ustar;
    //right shock wave
    double shock_speed;
    double rho_star2;
    //left rarefaction wave
    double c_star1;
    double rho_star1;
    double rare_speed_head;
    double rare_speed_tail;

public:
    Solve_System()
        :Pstar(0.0),Ustar(0.0),shock_speed(0.0),rho_star2(0.0)\
        ,c_star1(0.0),rho_star1(0.0),rare_speed_head(0.0),rare_speed_tail(0.0)
    {
    }
    double get_Pstar(){
        return Pstar;
    }
    double get_Ustar(){
        return Ustar;
    }
    double get_shock_speed(){
        return shock_speed;
    }
    double get_rho_star2(){
        return rho_star2;
    }
    double get_c_star1(){
        return c_star1;
    }
    double get_rho_star1(){
        return rho_star1;
    }
    double get_rare_speed_head(){
        return rare_speed_head;
    }
    double get_rare_speed_tail(){
        return rare_speed_tail;
    }

    void solve(const double Residual){
        int k = 0;
        double error = 0.0;
        double Pstar_old = Pstar;
        double Pstar_new = Pstar;
        double alpha_1_new = alpha_1_init;
        double alpha_2_new = alpha_2_init;

        k++;
        Pstar_new = (alpha_1_new*P_2 + alpha_2_new*P_1 - alpha_1_new*alpha_2_new*(u_1-u_2))\
                /(alpha_1_new + alpha_2_new);

        error = fabs(Pstar_new - Pstar_old);

        while(error > Residual){
            k++;
            Pstar_old = Pstar_new;
            alpha_1_new = rho_1*c_1*(1-Pstar_new/P_1)/((2*gamma/(gamma-1))*(1-std::pow(Pstar_new/P_1,(gamma-1)/(2*gamma))));
            alpha_2_new = rho_2*c_2*sqrt(((gamma+1)/(2*gamma))*(Pstar_new/P_2 -1)+1);
            Pstar_new = (alpha_1_new*P_2 + alpha_2_new*P_1 - alpha_1_new*alpha_2_new*(u_1-u_2))\
                    /(alpha_1_new + alpha_2_new);
            error = fabs(Pstar_new - Pstar_old);
            std::cout << "iters = " << k << " Pstar = " << Pstar_new << " error = " << error << std::endl;
        }
        Pstar = Pstar_new;
        Ustar = (alpha_1_new*u_1 + alpha_2_new*u_2 - (P_2-P_1))/(alpha_1_new + alpha_2_new);
        std::cout << "Pstar = " << Pstar << std::endl;
        std::cout << "Ustar = " << Ustar << std::endl;
        //right shock wave:
        shock_speed = u_2 + alpha_2_new/rho_2;
        rho_star2 = alpha_2_new/(shock_speed-Ustar);
        std::cout << "shock speed = " << shock_speed << std::endl;
        std::cout << "rho_star2 = " << rho_star2 << std::endl;

        //left rarefaction wave:
        c_star1 = c_1 + (gamma-1)*(u_1-Ustar)/2;
        rho_star1 = gamma*Pstar/(c_star1*c_star1);
        std::cout << "rho_star1 = " << rho_star1 << std::endl;
        rare_speed_head = u_1 - c_1;
        std::cout << "rarefaction head speed = " << rare_speed_head << std::endl;
        rare_speed_tail = Ustar - c_star1;
        std::cout << "rarefaction tail speed = " << rare_speed_tail << std::endl;
    }

};

int main(int, char**){
    //step 1: calculate the p_star and u_star
    Solve_System solver;
    solver.solve(1e-6);
    //step 2: calculate the correlated features
    // left is rarefaction wave
    // test accuracy
    /*
    double p_star = 0.0;
    p_star = solve.get_Pstar();
    double u_star = 0.0;
    u_star = solve.get_Ustar();
    double z_2 = u_2 - (p_star - P_2)/(rho_2*(u_2-u_star));
    double rho_star = rho_2*(u_2 - z_2)/(u_star - z_2);
    std::cout << "z_2 = " << z_2 << std::endl;
    std::cout << "rho_star = " << rho_star << std::endl;
    */
    //step 3: global solution
    const int N = 201;
    double u[N] = {0.0};
    double rho[N] = {0.0};
    double P[N] = {0.0};
    double time = 0.0;
    std::cout << "Please input the time you want to calculate: " << std::endl;
    std::cin >> time;
    if(time == 0.0){
        std::cout << "time can not be zero!" << std::endl;
        return 0;
    }
    for(int i=0;i<N;i++){
        double x = -1.0 + 2/(N-1)*i;
        double omega = x / time;
        if(omega <  solver.get_rare_speed_head()){

        }
    }





    return 0;
}
