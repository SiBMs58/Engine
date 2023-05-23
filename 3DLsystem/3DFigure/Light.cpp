//
// Created by Siebe Mees on 22/05/2023.
//

#include "Light.h"

Light::Light(const vector<double> &ambient, const vector<double> &diffuse, const vector<double> &specular) {
    ambientLight.red = ambient[0];
    ambientLight.green = ambient[1];
    ambientLight.blue = ambient[2];
    diffuseLight.red = diffuse[0];
    diffuseLight.green = diffuse[1];
    diffuseLight.blue = diffuse[2];
    specularLight.red = specular[0];
    specularLight.green = specular[1];
    specularLight.blue = specular[2];
}
