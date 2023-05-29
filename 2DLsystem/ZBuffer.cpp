//
// Created by Siebe Mees on 21/04/2023.
//

#include "ZBuffer.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>


ZBuffer::ZBuffer(const int width, const int height) {
    double posInf = std::numeric_limits<double>::infinity();
    for (int i = 0; i < width; ++i) {
        // Columns
        vector<double> col;
        for (int j = 0; j < height; ++j) {
            // Row
            col.push_back(posInf);
        }
        this->push_back(col);
    }
}

void ZBuffer::draw_zbuf_line(ZBuffer &zbuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, Color &color) {

    img::Color c = img::Color(color.red, color.green, color.blue);

    /**
     * a = length of the line
     * s = i, according to cursus
     * i = iterator of drawing the line (based on draw_line function)
     */

    if (x0 >= image.get_width() || y0 >= image.get_height() || x1 >= image.get_width() || y1 > image.get_height()) {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
           << image.get_width() << " and height " << image.get_height();
        throw std::runtime_error(ss.str());
    }
    if (x0 == x1)
    {
        //special case for x0 == x1
        double a = std::max(y0, y1) - std::min(y0, y1)+1;
        double s = a;
        unsigned int i = std::min(y0, y1);
        while (i <= std::max(y0, y1)) {
            double p = s / a;
            if (y0 > y1){
                std::swap(z0, z1);
                std::swap(x0, x1);
            }
            double Zi = p / z0 + (1 - p) / z1;
            if (Zi < zbuffer[x0][i]) {
                image(x0, i) = c;
                zbuffer[x0][i] = Zi;
            }
            i++;
            s--;
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        double a = std::max(x0, x1) - std::min(x0, x1)+1;
        double s = a;
        unsigned int i = std::min(x0, x1);
        while (i <= std::max(x0, x1)) {
            double p = s/a;
            if (x0 > x1) {
                std::swap(z0, z1);
                std::swap(y0, y1);
            }
            double Zi = p/z0 + (1-p)/z1;
            if (Zi < zbuffer[i][y0]) {
                image(i, y0) = c;
                zbuffer[i][y0] = Zi;
            }
            i++;
            s--;
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            double a = (x1 - x0)+1;
            double s = a;
            unsigned int i = 0;
            while (i <= (x1 - x0)) {
                double p = s/a;
                double Zi = p/z0 + (1-p)/z1;
                if (Zi < zbuffer[x0 + i][(unsigned int) round(y0 + m * i)]) {
                    image(x0 + i, (unsigned int) round(y0 + m * i)) = c;
                    zbuffer[x0 + i][(unsigned int) round(y0 + m * i)] = Zi;
                }
                s--;
                i++;
            }
        }
        else if (m > 1.0)
        {
            double a = (y1 - y0)+1;
            double s = a;
            unsigned int i = 0;
            while (i <= (y1 - y0)) {
                double p = i/a;
                double Zi = p/z0 + (1-p)/z1;
                if (Zi < zbuffer[(unsigned int) round(x0 + (i / m))][y0 + i]) {
                    image((unsigned int) round(x0 + (i / m)), y0 + i) = c;
                    zbuffer[(unsigned int) round(x0 + (i / m))][y0 + i] = Zi;
                }
                i++;
                s--;
            }
        }
        else if (m < -1.0)
        {
            double a = (y0 - y1)+1;
            double s = a;
            unsigned int i = 0;
            while (i <= (y0 - y1)) {
                double p = i/a;
                double Zi = p/z0 + (1-p)/z1;
                if (Zi < zbuffer[(unsigned int) round(x0 - (i / m))][y0 - i]) {
                    image((unsigned int) round(x0 - (i / m)), y0 - i) = c;
                    zbuffer[(unsigned int) round(x0 - (i / m))][y0 - i] = Zi;
                }
                i++;
                s--;
            }
        }
    }

}

void calculateXvals(const Point2D &p, const Point2D &q, const double y, double &xL, double &xR) {
    if ((y - p.y)*(y - q.y) <= 0 && p.y != q.y) {
        double xI = q.x + (p.x - q.x)*((y - q.y) / (p.y - q.y));
        if (xI < xL) {
            xL = xI;
        }
        if (xI > xR) {
            xR = xI;
        }
    }

}

Matrix eyePointTrans(Vector3D &eyepoint){
    double r = 0;
    double theta = 0;
    double phi =0;
    // 1. Omzetten van carthesische naar poolcoördinaten
    r = sqrt(pow(eyepoint.x, 2) + pow(eyepoint.y, 2) + pow(eyepoint.z, 2));
    theta = atan2(eyepoint.y, eyepoint.x);
    phi = acos((eyepoint.z/r));
    // 2. eyePointTransformationMatrix opstellen:
    Matrix eyePointTransformationMatrix;
    eyePointTransformationMatrix(1, 1) = -sin(theta);
    eyePointTransformationMatrix(1, 2) = -cos(theta)*cos(phi);
    eyePointTransformationMatrix(1, 3) = cos(theta)*sin(phi);
    eyePointTransformationMatrix(2, 1) = cos(theta);
    eyePointTransformationMatrix(2, 2) = -sin(theta)*cos(phi);
    eyePointTransformationMatrix(2, 3) = sin(theta)*sin(phi);
    eyePointTransformationMatrix(3,2) = sin(phi);
    eyePointTransformationMatrix(3,3) = cos(phi);
    eyePointTransformationMatrix(4,3) = -r;
    eyePointTransformationMatrix(4,4) = 1;
    return eyePointTransformationMatrix;
}

void diffuseColor(img::Color &c, int x, int y, const double &d, const double &dx, const double &dy, ZBuffer& zBuffer,  Lights3D &lights3D, Vector3D &normaalVector, const Color &diffuseReflection, const Color &specularReflection, const double &reflectionCoeff, Vector3D &eye) {

    double xproj = x; // - dx;
    double yproj = y; // - dy;

    double zEye = 1.0 / zBuffer[x][y];

    double xEye = -zEye * (xproj - dx) / d;
    double yEye = -zEye * (yproj - dy) / d;

    Vector3D point = Vector3D::point(xEye, yEye, zEye);

    Matrix eyePointTransformationMatrix = eyePointTrans(eye);
    Matrix invEyepointMatrix = Matrix::inv(eyePointTransformationMatrix);
    point *= invEyepointMatrix;

    for (Light *light: lights3D) {
        if (light->getType() == "PointLight") {

            Vector3D p = point * eyePointTransformationMatrix;
            Vector3D l = Vector3D::normalise(light->getLd() - p);
            double alpha = normaalVector.x * l.x + normaalVector.y * l.y + normaalVector.z * l.z;
            //double alpha = light->getAlpha(normaalVector);
            if (alpha > 0) {
                c.red = c.red + diffuseReflection.red * light->diffuseLight.red * alpha;
                c.green = c.green + diffuseReflection.green * light->diffuseLight.green * alpha;
                c.blue = c.blue + diffuseReflection.blue * light->diffuseLight.blue * alpha;
            }


            /*if (reflection) {

                Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - l);
//                Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - eLightPos);
                Vector3D pointToEye = Vector3D::normalise(eyepoint * eyePointM - ePoint);

                double cosBeta = r.x * pointToEye.x + r.y * pointToEye.y + r.z * pointToEye.z;

//                std::cout << cosBeta << std::endl;

                if (cosBeta > 0) {
                    color.red = saturate(color.red + specularReflection.red * light->specularLight.red *
                                                     pow(cosBeta, reflectionCoeff));
                    color.green = saturate(color.green + specularReflection.green * light->specularLight.green *
                                                         pow(cosBeta, reflectionCoeff));
                    color.blue = saturate(color.blue + specularReflection.blue * light->specularLight.blue *
                                                       pow(cosBeta, reflectionCoeff));
                }
            }*/
        }
        /*if (light->getType() == "INFLIGHT" and reflection) {

            // Alles in oorspronkelijke coordinaten
            Vector3D lightPos = light->getLd();
            Vector3D eLightPos = Vector3D::normalise(lightPos * TransMatrix::invEyePointTrans(eyepoint));
            Vector3D ePoint = point; // * TransMatrix::eyePointTrans(eyepoint);

//            Vector3D l = Vector3D::normalise(eLightPos-ePoint); // de light->getLd() en point zijn allebei de orginele posities
            Vector3D normaalVectorB = Vector3D::normalise(normaalVector * TransMatrix::invEyePointTrans(eyepoint));
            double alpha =
                    normaalVectorB.x * eLightPos.x + normaalVectorB.y * eLightPos.y + normaalVectorB.z * eLightPos.z;


//            Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - l);
            Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - eLightPos);
            Vector3D pointToEye = Vector3D::normalise(eyepoint - ePoint);

            double cosBeta = r.x * pointToEye.x + r.y * pointToEye.y + r.z * pointToEye.z;

//                std::cout << cosBeta << std::endl;

            if (cosBeta > 0) {
                color.red = saturate(color.red + specularReflection.red * light->specularLight.red *
                                                 pow(cosBeta, reflectionCoeff));
                color.green = saturate(color.green + specularReflection.green * light->specularLight.green *
                                                     pow(cosBeta, reflectionCoeff));
                color.blue = saturate(color.blue + specularReflection.blue * light->specularLight.blue *
                                                   pow(cosBeta, reflectionCoeff));
            }
        }*/
    }


}

void ZBuffer::draw_zbuf_triag(ZBuffer &zbuffer, img::EasyImage &image, const Vector3D &A, const Vector3D &B, const Vector3D &C,
                              double d, double dx, double dy, Color ambientReflection, Color diffuseReflection, Color specularReflection, double reflectionCoefficient, Lights3D& lights, Vector3D& eye) {

    // 1. Projectie van de driehoek
    Point2D projectA;
    double xA = (d*A.x/-A.z)+dx;
    double yA = (d*A.y/-A.z)+dy;
    projectA.x = xA;
    projectA.y = yA;
    Point2D projectB;
    double xB = (d*B.x/-B.z)+dx;
    double yB = (d*B.y/-B.z)+dy;
    projectB.x = xB;
    projectB.y = yB;
    Point2D projectC;
    double xC = (d*C.x/-C.z)+dx;
    double yC = (d*C.y/-C.z)+dy;
    projectC.x = xC;
    projectC.y = yC;

    /*// Berekenen van de normaalvector
    Vector3D AB = Vector3D::normalise(B - A);
    Vector3D AC = Vector3D::normalise(C - A);
    Vector3D n = Vector3D::cross(AB, AC);
    n = Vector3D::normalise(n);*/


    // 3.1 Coördinaten van het zwaartepunt van A’B’C’
    double xG = (projectA.x + projectB.x + projectC.x) / 3.0;
    double yG = (projectA.y + projectB.y + projectC.y) / 3.0;
    double zG = (1.0 / (3 * A.z)) + (1.0 / (3 * B.z)) + (1.0 / (3 * C.z));
    // 3.2 Berekening van dzdx, dzdy
    Vector3D u = B - A;
    Vector3D v = C - A;
    Vector3D w = Vector3D::cross(u, v);

    Color color;
    color.red = 0;
    color.green = 0;
    color.blue = 0;

    if (lights.empty()) {
        color.red = ambientReflection.red;
        color.green = ambientReflection.green;
        color.blue = ambientReflection.blue;
    }

    // Ambient licht wordt altijd toegepast
    for (Light *light: lights) {
        color.red += ambientReflection.red * light->ambientLight.red;
        color.green += ambientReflection.green * light->ambientLight.green;
        color.blue += ambientReflection.blue * light->ambientLight.blue;
    }
    for (auto *light: lights) {
        if (light->getType() == "InfLight") {
            // First applying the diffuse light
            Vector3D lOld = light->getLd(); // Creating l vector
            lOld.x = -lOld.x; lOld.y = -lOld.y; lOld.z = -lOld.z; // Creating l with the right values
            Vector3D n = Vector3D::normalise(w); // Normalising
            Vector3D l = Vector3D::normalise(lOld);
            //double divValue = sqrt(pow(w.x,2) + pow(w.y,2) + pow(w.z,2)); // Manually normalisation
            //n.x = n.x / divValue; n.y = n.y / divValue; n.z = n.z / divValue; // Manual normalisation
            double cosAlpha = n.x * l.x + n.y * l.y + n.z * l.z; // Calculating cosAlpha
            if (cosAlpha > 0) {
                // Else no diffuse light :(, but if cosAlpha > 0 we :), applying diffuse light
                color.red += light->diffuseLight.red * diffuseReflection.red * cosAlpha;
                color.green += light->diffuseLight.green * diffuseReflection.green * cosAlpha;
                color.blue += light->diffuseLight.blue * diffuseReflection.blue * cosAlpha;
            }

            // Now applying specular light
            // Defining vector r
            Vector3D r = Vector3D::vector(2 * cosAlpha * n.x - l.x,
                                          2 * cosAlpha * n.y - l.y,
                                          2 * cosAlpha * n.z - l.z);


        }  // If the light is not a point or infinity, leave it out, this light came into this engine ILLEGALLY
    }

    double k = (w.x * A.x) + (w.y * A.y) + (w.z * A.z);
    double dzdx = w.x / (-d * k);
    double dzdy = w.y / (-d * k);

    // 2. Bepalen van de pixels die tot de 3-hoek A’B’C’ behoren
    // 2.1 y-waardes
    for (int i = round(std::min({projectA.y, projectB.y, projectC.y})+0.5); i <= round(std::max({projectA.y, projectB.y, projectC.y})-0.5); i++) {
        double xL = std::numeric_limits<double>::infinity();
        double xR = -std::numeric_limits<double>::infinity();
        calculateXvals(projectA, projectB, i, xL, xR);
        calculateXvals(projectA, projectC, i, xL, xR);
        calculateXvals(projectB, projectC, i, xL, xR);
        /*// 3. De 1/z waarde voor elke pixel berekenen
        // 3.1 Coördinaten van het zwaartepunt van A’B’C’
        double xG = (projectA.x + projectB.x + projectC.x) / 3.0;
        double yG = (projectA.y + projectB.y + projectC.y) / 3.0;
        double zG = (1.0 / (3 * A.z)) + (1.0 / (3 * B.z)) + (1.0 / (3 * C.z));
        // 3.2 Berekening van dzdx, dzdy
        Vector3D u = B - A;
        Vector3D v = C - A;
        Vector3D w = Vector3D::cross(u, v);
        double k = (w.x * A.x) + (w.y * A.y) + (w.z * A.z);
        double dzdx = w.x / (-d * k);
        double dzdy = w.y / (-d * k);*/
        // 3.3 Voor alle pixel (x,y) element van 3-hoek A’B’C’:
        for (int j = xL + 0.5; j <= xR - 0.5; j++) {
            double z = (1.0001 * (zG)) + ((j - xG) * dzdx) + ((i - yG) * dzdy);
            if (z < zbuffer[j][i]) {
                Color saveColor = color;
                zbuffer[j][i] = z;
                for (auto *light: lights) {
                    if (light->getType() == "PointLight") {
                        double lX = ((j - dx) * -(1.0 / z)) / d;
                        double lY = ((i - dy) * -(1 / z)) / d;
                        double lZ = 1.0 / z;
                        Vector3D ld = Vector3D::vector(light->getLd().x, light->getLd().y, light->getLd().z);
                        ld.x -= lX;
                        ld.y -= lY;
                        ld.z -= lZ;
                        ld = Vector3D::normalise(ld);
                        Vector3D n = Vector3D::normalise(w);
                        double cosAlpha = ld.x * n.x + ld.y * n.y + ld.z * n.z;
                        if (cosAlpha > 0 and cosAlpha > cos(light->getSpotAngle())) {
                            double final = (cosAlpha - cos(light->getSpotAngle()))/(1-cos(light->getSpotAngle()));
                            color.red += light->diffuseLight.red * diffuseReflection.red * final;
                            color.green += light->diffuseLight.green * diffuseReflection.green * final;
                            color.blue += light->diffuseLight.blue * diffuseReflection.blue * final;
                        }
                        Vector3D rOld = Vector3D::vector(2 * cosAlpha * n.x - ld.x,
                                                         2 * cosAlpha * n.y - ld.y,
                                                         2 * cosAlpha * n.z - ld.z);
                        Vector3D r = Vector3D::normalise(rOld);
                        Vector3D m = Vector3D::vector(-lX,-lY,-lZ);
                        m = Vector3D::normalise(m);
                        double cosBeta = r.x * m.x + r.y * m.y + r.z * m.z;
                        // Doing the cosine of beta to the power of M_s
                        double cosBetaMS = pow(cosBeta, reflectionCoefficient);
                        if (cosBeta > 0) {
                            color.red += light->specularLight.red * specularReflection.red * cosBetaMS;
                            color.green += light->specularLight.green * specularReflection.green * cosBetaMS;
                            color.blue += light->specularLight.blue * specularReflection.blue * cosBetaMS;
                        } if (color.red > 255) {
                            color.red = 255;
                        } if (color.green > 255) {
                            color.green = 255;
                        } if (color.blue > 255) {
                            color.blue = 255;
                        }
                    }
                }
                img::Color finalColor = img::Color(color.red, color.green, color.blue);
                image(j, i) = finalColor;
                color = saveColor;
            }
        }
    }
}