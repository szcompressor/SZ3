//
// Created by Kai Zhao on 9/1/20.
//

#ifndef SZ_INTERPOLATORS_HPP
#define SZ_INTERPOLATORS_HPP

namespace SZ3 {
template <class T>
inline T interp_linear(T a, T b) {
    return (a + b) / 2;
}

template <class T>
inline T interp_linear1(T a, T b) {
    return -0.5 * a + 1.5 * b;
}

template <class T>
inline T interp_quad_1(T a, T b, T c) {
    return (3 * a + 6 * b - c) / 8;
}

template <class T>
inline T interp_quad_2(T a, T b, T c) {
    return (-a + 6 * b + 3 * c) / 8;
}

template <class T>
inline T interp_quad_3(T a, T b, T c) {
    return (3 * a - 10 * b + 15 * c) / 8;
}

template <class T>
inline T interp_cubic(T a, T b, T c, T d) {
    return (-a + 9 * b + 9 * c - d) / 16;
}

template <class T>
inline T interp_cubic_front(T a, T b, T c, T d) {
    return (5 * a + 15 * b - 5 * c + d) / 16;
}

template <class T>
inline T interp_cubic_front_2(T a, T b, T c, T d) {
    return (a + 6 * b - 4 * c + d) / 4;
}

template <class T>
inline T interp_cubic_back_1(T a, T b, T c, T d) {
    return (a - 5 * b + 15 * c + 5 * d) / 16;
}

template <class T>
inline T interp_cubic_back_2(T a, T b, T c, T d) {
    return (-5 * a + 21 * b - 35 * c + 35 * d) / 16;
}

template <class T>
inline T interp_cubic2(T a, T b, T c, T d) {
    return (-3 * a + 23 * b + 23 * c - 3 * d) / 40;
}

template <class T>
inline T interp_akima(T a, T b, T c, T d) {
    T t0 = 2 * b - a - c;
    T t1 = 2 * c - b - d;
    T abt0 = fabs(t0);
    T abt1 = fabs(t1);
    if (fabs(abt0 + abt1) > 1e-9) {
        return (b + c) / 2 + (t0 * abt1 + t1 * abt0) / 8 / (abt0 + abt1);
    } else {
        return (b + c) / 2;
    }
}

template <class T>
inline T interp_pchip(T a, T b, T c, T d) {
    T pchip = (b + c) / 2;
    if ((b - a < 0) == (c - b < 0) && fabs(c - a) > 1e-9) {
        pchip += 1 / 4 * (b - a) * (c - b) / (c - a);
    }
    if ((c - b < 0) == (d - c < 0) && fabs(d - b) > 1e-9) {
        pchip -= 1 / 4 * (c - b) * (d - c) / (d - b);
    }
    return pchip;
}
    namespace QoZ {
        template<class T>
        inline T interp_linear(T a, T b) {
            return (a + b) / 2;
        }

        template<class T>
        inline T interp_linear1(T a, T b) {
            return -0.5 * a + 1.5 * b;
        }

        template<class T>
        inline T interp_quad_1(T a, T b, T c) {
            return (3 * a + 6 * b - c) / 8;
        }
        template<class T>
        inline T interp_quad_1_adj(T a, T b, T c) {
            return ( a + 3 * b - c) / 3;
        }

        template<class T>
        inline T interp_quad_2(T a, T b, T c) {
            return (-a + 6 * b + 3 * c) / 8;
        }

        template<class T>
        inline T interp_quad_2_adj(T a, T b, T c) {
            return (-a + 3 * b +  c) / 3;
        }

        template<class T>
        inline T interp_quad_3(T a, T b, T c) {
            return (3 * a - 10 * b + 15 * c) / 8;
        }

        template<class T>
        inline T interp_quad_3_adj(T a, T b, T c) {
            return 3 * c - 3 * b + a;
        }

        template<class T>
        inline T interp_cubic_1(T a, T b, T c, T d) {
            return (-a + 9 * b + 9 * c - d) / 16;//noknot
            //return -0.06368435202786181*a+0.5731591682507563*b+0.5731591682507563*c-0.06368435202786181*d;
            //return (-3*a+23*b+23*c-3*d)/40;
        }

        template<class T>
        inline T interp_cubic_2(T a, T b, T c, T d) {
            return (-3 * a + 23 * b + 23 * c - 3 * d) / 40;//nat
        }

        template<class T>
        inline T interp_cubic(uint8_t cst, T a, T b, T c, T d){
            if (cst==0)
                return (-a + 9 * b + 9 * c - d) / 16;
            else
                return (-3 * a + 23 * b + 23 * c - 3 * d) / 40;
        }

       

        template<class T>
        inline T interp_cubic_adj(uint8_t cst, T a, T b, T c, T d,T e,T f) {
            if (cst==0)
                return (-b+4*c+4*d-e)/6;
            else
                return (3*a-18*b+46*c+46*d-18*e+3*f)/62;
        }

        template<class T>
        inline T interp_cubic_adj2(uint8_t cst, T a, T b, T c, T d,T f) {
            if (cst==0)
                return (-4*b+15*c+10*d-f)/20;
            else
                return (12*a-72*b+181*c+118*d-15*f)/224;
        }

        





        template<class T>
        inline T interp_cubic_adj_1(T a, T b, T c, T d,T e,T f) {//adj6 nat
            return (3*a-18*b+46*c+46*d-18*e+3*f)/62;
            //return (-3*b+11*c+11*d-3*e)/16;
        }
        template<class T>
        inline T interp_cubic_adj_2(T a, T b, T c, T d,T e,T f) {//adj6 noknot
            return (-b+4*c+4*d-e)/6;
        }
        template<class T>
        inline T interp_cubic_adj_3(T a, T b, T c, T d,T f) {//adj5 nat
            return (12*a-72*b+181*c+118*d-15*f)/224;
            //return (-3*b+11*c+11*d-3*e)/16;
        }
        template<class T>
        inline T interp_cubic_adj_4(T a, T b, T c, T d,T f) {//adj5 noknot
            return (-4*b+15*c+10*d-f)/20;
            //return (-3*b+11*c+11*d-3*e)/16;
        }


        template<class T>
        inline T interp_2d(T a, T b, T c, T d) {
            return ( a + b + c + d) / 4;
        }

        template<class T>
        inline T interp_3d(T a, T b, T c, T d, T e,T f) {
            return ( a + b + c + d + e + f ) / 6;
        }

        template<class T>
        inline T lorenzo_1d(T a, T b) {
            return 2*b-a;
        }

        template<class T>
        inline T lorenzo_2d(T a, T b, T c) {
            return (b+c-a);
        }
        
        template<class T>
        inline T lorenzo_3d(T a, T b, T c, T d, T e,T f,T g) {
            return (a-b-c+d-e+f+g);
        }

        template<class T>
        inline T interp_ave3(T a, T b, T c) {
            return (a + b+c) / 3;
        }

        template<class T>
        inline T interp_cubic_front(T a, T b, T c, T d) {
            return (5 * a + 15 * b - 5 * c + d) / 16;
        }

        template<class T>
        inline T interp_cubic_front_adj(T a, T b, T c, T d) {
            return (17 * a + 44 * b - 18 * c + 3*d) / 46;
        }

        template<class T>
        inline T interp_cubic_front_2(T a, T b, T c, T d) {
            return ( a + 6 * b - 4 * c + d) / 4;
        }

        template<class T>
        inline T interp_cubic_back_1(T a, T b, T c, T d) {
            return (a - 5 * b + 15 * c + 5 * d) / 16;
        }

        template<class T>
        inline T interp_cubic_back_adj(T a, T b, T c, T d) {
            return (3*a - 18 * b + 44 * c + 17 * d) / 46;
        }

        template<class T>
        inline T interp_cubic_back_2(T a, T b, T c, T d) {
            return (-5 * a + 21 * b - 35 * c + 35 * d) / 16;
        }

        

        /*
        template<class T>
        inline T lanczos(T x, int a) {
            if(x==0)
                return 0;
            else if (fabs(x)>a)
                return 1;
            else{
                T pix=M_PI*x;
                return a*sin(pix)*sin(pix/a)/(pix*pix);
            }
        }
        template<class T>
        inline T interp_lanczos_2(T a, T b, T c, T d) {

            return a*lanczos(1.5,2)+b*lanczos(0.5,2)+c*lanczos(-0.5,2)+d*lanczos(-1.5,2);
            
            -0.06368435202786181
            0.5731591682507563

            
        }
        */

    }



}  // namespace SZ3
#endif  // SZ_INTERPOLATORS_HPP
