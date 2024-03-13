#include "bmp388.hpp"
#include <cstdio>
#include "pico/stdlib.h"
#include "hardware/i2c.h"
#include <stdio.h>
#include "pico/time.h"
#include <cmath>
#include "tusb.h"


#define I2C_PORT i2c0
#define I2C_SDA 12
#define I2C_SCL 13

BMP388 alt(I2C_PORT);



int main() {
    stdio_init_all();

    i2c_init(I2C_PORT, 400 * 1000);

    gpio_set_function(I2C_SDA, GPIO_FUNC_I2C);
    gpio_set_function(I2C_SCL, GPIO_FUNC_I2C);

    gpio_pull_up(I2C_SDA);
    gpio_pull_up(I2C_SCL);

#ifdef VERBOSE
    while (!tud_cdc_connected()) {
        sleep_ms(500);
    }
    printf("Connected\n");
#endif

    bool begun = alt.begin();
    double altitude; // variable actually containing altitude
    bool ret; // boolean checking to see if reading is outputted
    int points_read = 0; // you have read zero altitude points so far
    int alt_points;
    int vel_points;
    int SMA_vel_points;
    int total = 0;
    int SMA_length = 256;
    double prediction;
    int prev_altitude;
    int velocity;

    // variables for density calculations
    int mass = 120; // mass of LV dry in pounds
    double Cd = 0.53; // drag coefficient of LV
    double rho; // air density
    double A = 3.141592*3*3; // cross-sectional area of LV
    double a = -0.0065; // constant in density equation
    double g = 9.81; // gravitational constant
    double Tref = 288.15; // temperature in Kelvin
    double rhoref = 1.225; // reference density
    double R = 287; // ideal gas constant

    // varibales for euler's method
    double eul_pos;
    double eul_vel;
    double eul_acc;
    double timestep = 0.01;

    ret = alt.read_altitude(&altitude, 1013.25);
    while (altitude < 1000){
        ret = alt.read_altitude(&altitude, 1013.25);
        printf("alt: %.3f\n", altitude);
    } // will exit altitude loop when crosses 2000 feet


    while (true) {

        for (int i = 0; i < SMA_length+2; i++) {
            ret = alt.read_altitude(&altitude, 1013.25); // checking to see if altimeter is outputting a reading
    
            if (!ret) {
                printf("reading not successful\n"); 

            } else {

                alt_points[&points_read] = altitude; // plug in current altitude into altitude vector
                printf("alt: %.3f\n", altitude);

                if (points_read > 0) { // obtain first velocity reading b/c >1 altimeter measurements now exist
                    vel_points[&i - 1] = (altitude - alt_points[&i-1]); // need to find time between points idk how
                    total = total + vel_points[&i - 1];
                }
            }
        }
        // loop is exited after 257 points
        points_read = SMA_length+1; // 257 points were read in that for loop, so 256 velocity points are currently in the velocity vector.

        int j = 0;
        double SMA = total / SMA_length; // obtain velocity moving average
        SMA_vel_points[&j] = 0.9765*SMA - 27.1192; // adjust lagging velocity SMA to truer value. Must be changed for 256 adjustment
        // might not be necessary to log SMA into a vector
        SMA = 0.9765*SMA - 27.1192; // potentially just use this

        // maybe insert an altitude check so it doesn't enter into this phase until altitude meets height.

        int lru = 0;

        while (true) { // enter into predictive phase after first velocity SMA is obtained
            // dynamics phase
            eul_pos = altitude; // initial condition altitude
            eul_vel = SMA; // initial condition velocity
            rho = rhoref*pow((1+(a*eul_pos)/Tref),-g/(a*R)-1);
            

            while (eul_vel > 0) { // you've hit max altitude when euler velocity goes negative
                eul_acc = -g - 0.5*Cd*rho*eul_vel*eul_vel*A/mass; // obtain acceleration
                eul_vel = eul_vel + eul_acc*timestep; // update velocity 
                eul_pos = eul_pos + eul_pos*timestep; // update position
                rho = rhoref*pow((1+(a*eul_pos)/Tref),-g/(a*R)-1); // update density
            }

            prediction = eul_pos; // prediction is position when velocity goes negative (apogee)

            // Log altitude and prediction to SD card here (literally no idea how to do this rn) //
            
            prev_altitude = altitude;
            // Updating measurements and SMA
            ret = alt.read_altitude(&altitude, 1013.25); // obtain new altitude

            velocity = (altitude - prev_altitude);  // obtain new velocity
            // figure out time between measurements somehow

            total = total + velocity - vel_points[&lru]; // add newest measurement to total and subtract oldest
            SMA = total / SMA_length; // update SMA

            // Must figure out how to shift all velocity measurements to the left one spot in lines below
            vel_points[&lru] = velocity; // replace newest velocity
            lru = (lru + 1) % SMA_length; // LRU stuff

        }
    }
    return 0;
    
}