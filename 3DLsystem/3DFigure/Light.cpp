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

string Light::getType() {
    return "Light";
}

Vector3D Light::getLd() {
    return Vector3D();
}

double Light::getAlpha(Vector3D &ld) {
    return 0;
}

double Light::getSpotAngle() {
    return 0;
}

InfLight::InfLight(const vector<double> &ambient, const vector<double> &diffuse, const vector<double> &specular, const Vector3D &ld) : Light(ambient, diffuse, specular) {
    ldVector = ld;
}

Vector3D InfLight::getLd() {
    return ldVector;
}

string InfLight::getType() {
    return "InfLight";
}

double InfLight::getAlpha(Vector3D &ld) {
    return ldVector.x * ld.x + ldVector.y * ld.y + ldVector.z * ld.z;
}

double InfLight::getSpotAngle() {
    return 0;
}

PointLight::PointLight(const vector<double> &ambient, const vector<double> &diffuse, const vector<double> &specular, const Vector3D &loc) : Light(ambient, diffuse, specular){
    location = loc;
}

Vector3D PointLight::getLd() {
    return location;
}

string PointLight::getType() {
    return "PointLight";
}

double PointLight::getSpotAngle() {
    return spotAngle;
}

void PointLight::setSpotAngle(double spotAngle) {
    this->spotAngle = spotAngle;
}