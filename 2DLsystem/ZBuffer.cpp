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

void ZBuffer::draw_zbuf_triag(ZBuffer &zbuffer, img::EasyImage &image, const Vector3D &A, const Vector3D &B, const Vector3D &C,
                              double d, double dx, double dy, Color ambientReflection, Color diffuseReflection, Color specularReflection, double reflectionCoefficient, Lights3D& lights) {

    img::Color c = img::Color(ambientReflection.red, ambientReflection.green, ambientReflection.blue);

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

    // 2. Bepalen van de pixels die tot de 3-hoek A’B’C’ behoren
    // 2.1 y-waardes
    for (int i = round(std::min({projectA.y, projectB.y, projectC.y})+0.5); i <= round(std::max({projectA.y, projectB.y, projectC.y})-0.5); i++) {
        double xL = std::numeric_limits<double>::infinity();
        double xR = -std::numeric_limits<double>::infinity();
        calculateXvals(projectA, projectB, i, xL, xR);
        calculateXvals(projectA, projectC, i, xL, xR);
        calculateXvals(projectB, projectC, i, xL, xR);
        // 3. De 1/z waarde voor elke pixel berekenen
        // 3.1 Coördinaten van het zwaartepunt van A’B’C’
        double xG = (projectA.x+projectB.x+projectC.x)/3.0;
        double yG = (projectA.y+projectB.y+projectC.y)/3.0;
        double zG = (1.0/(3*A.z))+(1.0/(3*B.z))+(1.0/(3*C.z));
        // 3.2 Berekening van dzdx, dzdy
        Vector3D u = B-A;
        Vector3D v = C-A;
        Vector3D w = Vector3D::cross(u, v);
        double k = (w.x*A.x)+(w.y*A.y)+(w.z*A.z);
        double dzdx = w.x/(-d*k);
        double dzdy = w.y/(-d*k);
        // 3.3 Voor alle pixel (x,y) element van 3-hoek A’B’C’:
        for (int j = xL+0.5; j <= xR-0.5; j++) {
            double z = (1.0001*(zG))+((j-xG)*dzdx)+((i-yG)*dzdy);
            if (z < zbuffer[j][i]) {
                image(j, i) = c;
                zbuffer[j][i] = z;
            }
        }
    }



}