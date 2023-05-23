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

};

class InfLight: public Light
{
public:
    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;
};

class PointLight: public Light
{
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;
};



typedef vector<Light> Lights3D;


#endif //ENGINE_LIGHT_H
