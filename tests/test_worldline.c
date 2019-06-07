#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

#include "../Thirring.h"
int **field;
int **diraclink;


/* Create configurations with different negative loops for testing */
void loop_around_mass_at(int t, int x){
    field[t+1][x+1] = MONOMER;
    diraclink[t+1][x+1] = NDIRS;
    diraclink[t+0][x+0] = TUP;
    diraclink[t+1][x+0] = TUP;
    diraclink[t+2][x+0] = XUP;
    diraclink[t+2][x+1] = XUP;
    diraclink[t+2][x+2] = TDN;
    diraclink[t+1][x+2] = TDN;
    diraclink[t+0][x+2] = XDN;
    diraclink[t+0][x+1] = XDN;
}

void mass_without_loop_at(int t, int x){
    field[t+1][x+1] = MONOMER;
    diraclink[t+1][x+1] = NDIRS;
    diraclink[t+0][x+0] = TUP;
    diraclink[t+1][x+0] = TDN;
    diraclink[t+2][x+0] = XUP;
    diraclink[t+2][x+1] = XDN;
    diraclink[t+2][x+2] = TDN;
    diraclink[t+1][x+2] = TUP;
    diraclink[t+0][x+2] = XDN;
    diraclink[t+0][x+1] = XUP;
}

void negative_time_loop(){
    /* Write over the whole config */
    for(int t=0; t<NT; t++) for(int x=0;x<NX;x++){
        diraclink[t][x]= XUP*(1-x%2) + XDN*(x%2);
    }
    /* Time wrapping loop */
    for(int t=0; t<NT; t++){
        diraclink[t][2]= TUP*(1-t%2) + TDN*(t%2);
        diraclink[t][3]= TUP;
    }
    /* Add a kink and shift links around it */
    diraclink[0][3]= XDN;
    diraclink[0][2]= TUP;
    diraclink[1][2]= TUP;
    diraclink[2][2]= XUP;
    diraclink[1][0]= XDN;
    for(int x=3; x<NX; x++){
        diraclink[1][x]= XUP*(x%2) + XDN*(1-x%2);
    }
    diraclink[3][2]= XDN;
    diraclink[3][1]= XUP;
    diraclink[3][0]= TDN;
    diraclink[2][0]= TUP;
    diraclink[2][1]= TDN;
    diraclink[1][1]= TUP;
}

/* A test that checks if 1 is true */
static void test_sector(void **state) {
    /* Allocation and Setup */
    setup_lattice(12345);

    /* Simple sector 0 configuration */
    for(int t=0; t<NT; t++) for(int x=0;x<NX;x++){
        field[t][x]=0;
        diraclink[t][x]= XUP*(1-x%2) + XDN*(x%2);
    }

    int sector = count_negative_loops();
    int sign = configuration_sign();

    assert_true(sector == 0);
    assert_true(sign == 1-(sector%2)*2);

    /* A time wrapping loop */
    for(int t=0; t<NT; t++){
        diraclink[t][0]= TUP*(1-t%2) + TDN*(t%2);
        diraclink[t][1]= TUP;
    }
    assert_true(count_negative_loops() == 0);

    /* Add a negative corner */
    if( NX > 5 ){
        negative_time_loop();
        assert_true(count_negative_loops() == 1);
    }

    if( NX > 5 && NT > 2 ){
        /* Add two negative loop around a mass monomer */
        for(int t=0; t<NT; t++) for(int x=0;x<NX;x++){
            diraclink[t][x]= XUP*(1-x%2) + XDN*(x%2);
        }
        loop_around_mass_at(0,0);
        loop_around_mass_at(0,3);

        assert_true(count_negative_loops() == 2);

        /* Only loop around one of the masses */
        loop_around_mass_at(0,0);
        mass_without_loop_at(0,3);

        assert_true(count_negative_loops() == 1);
    }
}

/* In the main function create the list of the tests */
int main(void) {
   const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_sector),
   };

   // Call a library function that will run the tests
   return cmocka_run_group_tests(tests, NULL, NULL);
}
