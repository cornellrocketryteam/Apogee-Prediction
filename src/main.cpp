#include "bmp388.hpp"
#include <cstdio>
#include "pico/stdlib.h"
#include "hardware/i2c.h"


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

    bool begun = alt.begin();
    double altitude;
    bool ret;

    while (true) {
        ret = alt.read_altitude(&altitude, 1013.25);
        if (!ret) {
            printf("reading not successful\n");
        } else {
            printf("alt: %.3f\n", altitude);
        }
        
    }
    return 0;
}