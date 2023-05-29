//
// Created by Siebe Mees on 22/05/2023.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H

#include <vector>
#include "../../2DLsystem/2DLine/Color.h"
#include "../../vector3d.h"
using namespace std;

class Light {
public:
    Light(const vector<double> &ambientLight, const vector<double> &diffuseLight, const vector<double> &specularLight);

    //de ambiente licht component
    Color ambientLight;
    //de diffuse licht component
    Color diffuseLight;
    //de diffuse licht component
    Color specularLight;

    virtual Vector3D getLd();
    virtual string getType();
    virtual double getAlpha(Vector3D &ld);
    virtual double getSpotAngle();

};

class InfLight: public Light {
public:
    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;

    InfLight(const vector<double> &ambient, const vector<double> &diffuse, const vector<double> &specular, const Vector3D &ld);

    virtual Vector3D getLd();
    virtual string getType();
    virtual double getAlpha(Vector3D &ld);
    virtual double getSpotAngle();

};

class PointLight: public Light
{
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;

    PointLight(const vector<double> &ambient, const vector<double> &diffuse, const vector<double> &specular, const Vector3D &loc);

    virtual Vector3D getLd();
    virtual string getType();
    virtual double getSpotAngle();
    void setSpotAngle(double spotAngle);
};



typedef vector<Light*> Lights3D;


#endif //ENGINE_LIGHT_H
