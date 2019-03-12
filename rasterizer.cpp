#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "hw2_types.h"
#include "hw2_math_ops.h"
#include "hw2_file_ops.h"
#include <iostream>

void swapVector(Vec3 &, Vec3 &);
void midpointInterpolation(Vec3, Vec3);

void triangleRasterization(Vec3, Vec3, Vec3);
double calculatef(Vec3, Vec3, double, double);

Camera cameras[100];
int numberOfCameras = 0;

Model models[1000];
int numberOfModels = 0;

Color colors[100000];
int numberOfColors = 0;

Translation translations[1000];
int numberOfTranslations = 0;

Rotation rotations[1000];
int numberOfRotations = 0;

Scaling scalings[1000];
int numberOfScalings = 0;

Vec3 vertices[100000];
int numberOfVertices = 0;

Color backgroundColor;

// backface culling setting, default disabled
int backfaceCullingSetting = 0;

Color **image;

/*
    Initializes image with background color
*/
void initializeImage(Camera cam)
{
    int i, j;

    for (i = 0; i < cam.sizeX; i++)
        for (j = 0; j < cam.sizeY; j++)
        {
            image[i][j].r = backgroundColor.r;
            image[i][j].g = backgroundColor.g;
            image[i][j].b = backgroundColor.b;
        }
}

/*
    Transformations, culling, rasterization are done here.
    You can define helper functions inside this file (rasterizer.cpp) only.
    Using types in "hw2_types.h" and functions in "hw2_math_ops.cpp" will speed you up while working.
*/
void forwardRenderingPipeline(Camera cam)
{
    // TODO: IMPLEMENT HERE
    for (int i = 0; i < numberOfModels; i++)
    { // Apply for every model
        Model currentModel = models[i];
        int numTriangles = currentModel.numberOfTriangles; // Number of Triangles

        for (int j = 0; j < numTriangles; j++)
        { // Apply for every triangle
            Triangle currentTriangle = currentModel.triangles[j];
            //
            Vec3 triangleVertices[3];

            for (int k = 0; k < 3; k++)
            {                                                               // Apply for every vertex
                int vertexID = currentTriangle.vertexIds[k];                // Find vertexid
                Vec3 vertex = vertices[vertexID];                           // Take vertex position
                double vertexVector[4] = {vertex.x, vertex.y, vertex.z, 1}; // Vertex 4 member vector to use in transformations
                int numTrans = currentModel.numberOfTransformations;
                double vertexResult[4] = {0, 0, 0, 0};
                for (int l = 0; l < numTrans; l++)
                { //   Apply transformations
                    int transID = currentModel.transformationIDs[l];
                    char a = currentModel.transformationTypes[l];
                    if (a == 's')
                    { // Transformation is a scaling
                        // TODO: do scaling
                        Scaling currentScaling = scalings[transID];
                        double scaleMatrix[4][4] = {
                            {currentScaling.sx, 0, 0, 0},
                            {0, currentScaling.sy, 0, 0},
                            {0, 0, currentScaling.sz, 0},
                            {0, 0, 0, 1}};

                        multiplyMatrixWithVec4d(vertexResult, scaleMatrix, vertexVector);
                        for (int tt = 0; tt < 4; tt++)
                            vertexVector[tt] = vertexResult[tt]; // Vertex vector updated
                    }
                    else if (a == 't')
                    { // Transformation is a translation
                        // TODO: do translation
                        Translation currentTranslation = translations[transID];
                        double translationMatrix[4][4] = {
                            {1, 0, 0, currentTranslation.tx},
                            {0, 1, 0, currentTranslation.ty},
                            {0, 0, 1, currentTranslation.tz},
                            {0, 0, 0, 1}};

                        multiplyMatrixWithVec4d(vertexResult, translationMatrix, vertexVector);
                        for (int tt = 0; tt < 4; tt++)
                            vertexVector[tt] = vertexResult[tt]; // Vertex vector updated
                    }
                    else if (a == 'r')
                    { // Transformation is a rotation
                        // TODO: do rotation



                        Rotation currentRotation = rotations[transID];
                        Vec3 u;
                        u.x = currentRotation.ux;
                        u.y = currentRotation.uy;
                        u.z = currentRotation.uz;

            std::cout << "Rotation: u: " << u.x << "," << u.y << "," << u.z << std::endl;
                        Vec3 v;
                        // calculate v
                        double min = std::min(currentRotation.ux, std::min(currentRotation.uy, currentRotation.uz));
                        if (currentRotation.ux == min)
                        {
                            v.x = 0;
                            v.y = -1 * currentRotation.uz;
                            v.z = currentRotation.uy;
                        }
                        else if (currentRotation.uy == min)
                        {
                            v.y = 0;
                            v.x = -1 * currentRotation.uz;
                            v.z = currentRotation.ux;
                        }
                        else if (currentRotation.uz == min)
                        {
                            v.z = 0;
                            v.x = -1 * currentRotation.uy;
                            v.y = currentRotation.ux;
                        }
                        else
                            std::cout << "Error while creating v vector in rotation" << std::endl;
                        // calculate w
                        Vec3 w = crossProductVec3(u, v);
                        double M[4][4] = {
                            {u.x, u.y, u.z, 0},
                            {v.x, v.y, v.z, 0},
                            {w.x, w.y, w.z, 0},
                            {0, 0, 0, 1}};
                        double inverseM[4][4] = {
                            {u.x, v.x, w.x, 0},
                            {u.y, v.y, w.y, 0},
                            {u.z, v.z, w.z, 0},
                            {0, 0, 0, 1}};
                        double rotationMatrix[4][4] = {
                            {1, 0, 0, 0},
                            {0, cos(currentRotation.angle * M_PI / 180), -1 * sin(currentRotation.angle * M_PI / 180), 0},
                            {0, sin(currentRotation.angle * M_PI / 180), cos(currentRotation.angle * M_PI / 180), 0},
                            {0, 0, 0, 1}};
                        double result[4][4];
                        double rotationFinal[4][4];
                        multiplyMatrixWithMatrix(result, rotationMatrix, M);
                        multiplyMatrixWithMatrix(rotationFinal, inverseM, result);

                        multiplyMatrixWithVec4d(vertexResult, rotationFinal, vertexVector);
                        for (int tt = 0; tt < 4; tt++)
                            vertexVector[tt] = vertexResult[tt]; // Vertex vector updated

            std::cout << "vertexVector: " << vertexVector[0] << "," << vertexVector[1] << "," << vertexVector[2] << std::endl; 
                    }
                    else
                    {
                        std::cout << "transformation type is unknown" << std::endl;
                    }

                    //
                    //
                    //
                }
                //
                // Eye transformation
                //
                double camMatrix[4][4] = {
                    {cam.u.x, cam.u.y, cam.u.z, -1 * (cam.u.x * cam.pos.x + cam.u.y * cam.pos.y + cam.u.z * cam.pos.z)},
                    {cam.v.x, cam.v.y, cam.v.z, -1 * (cam.v.x * cam.pos.x + cam.v.y * cam.pos.y + cam.v.z * cam.pos.z)},
                    {cam.w.x, cam.w.y, cam.w.z, -1 * (cam.w.x * cam.pos.x + cam.w.y * cam.pos.y + cam.w.z * cam.pos.z)},
                    {0, 0, 0, 1}};

                multiplyMatrixWithVec4d(vertexResult, camMatrix, vertexVector);
                for (int tt = 0; tt < 4; tt++)
                    vertexVector[tt] = vertexResult[tt]; // Vertex vector updated
                //
                // Projection Transformation
                //
                double perspectiveMatrix[4][4] = {
                    {(2 * cam.n) / (cam.r - cam.l), 0, (cam.r + cam.l) / (cam.r - cam.l), 0},
                    {0, (2 * cam.n) / (cam.t - cam.b), (cam.t + cam.b) / (cam.t - cam.b), 0},
                    {0, 0, (cam.f + cam.n) / (cam.n - cam.f), (2 * cam.f * cam.n) / (cam.n - cam.f)},
                    {0, 0, -1, 0}};

                multiplyMatrixWithVec4d(vertexResult, perspectiveMatrix, vertexVector);
                for (int tt = 0; tt < 4; tt++)
                    vertexVector[tt] = vertexResult[tt]; // Vertex vector updated
                vertexVector[0] = vertexVector[0] / vertexVector[3];
                vertexVector[1] = vertexVector[1] / vertexVector[3];
                vertexVector[2] = vertexVector[2] / vertexVector[3];
                vertexVector[3] = vertexVector[3] / vertexVector[3];
                //
                // Viewport Transformation
                //
                double viewportMatrix[4][4] = {
                    {cam.sizeX * 0.5, 0, 0, (cam.sizeX - 1) * 0.5},
                    {0, cam.sizeY * 0.5, 0, (cam.sizeY - 1) * 0.5},
                    {0, 0, 0.5, 0.5},
                    {0, 0, 0, 0}};

                multiplyMatrixWithVec4d(vertexResult, viewportMatrix, vertexVector);
                for (int tt = 0; tt < 4; tt++)
                    vertexVector[tt] = vertexResult[tt]; // Vertex vector updated
                Vec3 temp;
                temp.x = vertexVector[0];
                temp.y = vertexVector[1];
                temp.z = vertexVector[2];
                temp.colorId = vertex.colorId;
                triangleVertices[k] = temp;
                Color color = colors[vertex.colorId];

                //
            }

            // Backface Culling

            Vec3 ba = subtractVec3(triangleVertices[1], triangleVertices[0]);
            Vec3 ca = subtractVec3(triangleVertices[2], triangleVertices[0]);
            Vec3 normal = normalizeVec3(crossProductVec3(ba, ca));
            Vec3 triCenter = multiplyVec3WithScalar(addVec3(triangleVertices[2], addVec3(triangleVertices[0], triangleVertices[1])), 0.333);
            Vec3 viewVector = normalizeVec3(subtractVec3(triCenter, cam.gaze));

            if (backfaceCullingSetting == 0 || dotProductVec3(normal, viewVector) >= 0)
            {

                // Mid Point Algorithm
                for (int mid = 0; mid < 3; mid++)
                {
                    Vec3 v0;
                    Vec3 v1;
                    if (mid == 0)
                    {
                        v0 = triangleVertices[0];
                        v1 = triangleVertices[1];
                    }
                    else if (mid == 1)
                    {
                        v0 = triangleVertices[1];
                        v1 = triangleVertices[2];
                    }
                    else if (mid == 2)
                    {
                        v0 = triangleVertices[2];
                        v1 = triangleVertices[0];
                    }

                    v0.x = (int)v0.x;
                    v0.y = (int)v0.y;
                    v0.z = (int)v0.z;
                    v1.x = (int)v1.x;
                    v1.y = (int)v1.y;
                    v1.z = (int)v1.z;

                    midpointInterpolation(v0, v1);
                }
                if (currentModel.type == 1)
                {
                    Vec3 v0 = triangleVertices[0];
                    Vec3 v1 = triangleVertices[1];
                    Vec3 v2 = triangleVertices[2];
                    triangleRasterization(v0, v1, v2);
                }
            }
            /*
            Color c;
            c.r = 0;
            c.g = 0;
            c.b = 0;

            Color c1;
            c1.r = 255;
            c1.g = 0;
            c1.b = 0;

            image[(int)triangleVertices[0].x][(int)triangleVertices[0].y] = c1;
            image[(int)triangleVertices[1].x][(int)triangleVertices[1].y] = c;
            image[(int)triangleVertices[2].x][(int)triangleVertices[2].y] = c;*/
        }
    }
}

void triangleRasterization(Vec3 v0, Vec3 v1, Vec3 v2)
{
    double x0 = v0.x;
    double x1 = v1.x;
    double x2 = v2.x;
    double y0 = v0.y;
    double y1 = v1.y;
    double y2 = v2.y;
    double maxy = std::max(v0.y, std::max(v1.y, v2.y));
    double maxx = std::max(v0.x, std::max(v1.x, v2.x));
    double miny = std::min(v0.y, std::min(v1.y, v2.y));
    double minx = std::min(v0.x, std::min(v1.x, v2.x));
    double alpha, beta, gamma;

    /*miny = floor(miny);
    minx = floor(minx);
    maxx = ceil(maxx);
    maxy = ceil(maxy);*/

    Color color;
    Color color0 = colors[v0.colorId];
    Color color1 = colors[v1.colorId];
    Color color2 = colors[v2.colorId];
    for (int y = miny; y <= maxy; y++)
    {
        for (int x = minx; x <= maxx; x++)
        {
            alpha = calculatef(v1, v2, x, y) / calculatef(v1, v2, x0, y0);
            beta = calculatef(v2, v0, x, y) / calculatef(v2, v0, x1, y1);
            gamma = calculatef(v0, v1, x, y) / calculatef(v0, v1, x2, y2);
            if ((alpha) >= 0 && (beta) >= 0 && (gamma) >= 0)
            {
                color.r = (alpha * color0.r + beta * color1.r + gamma * color2.r);
                color.g = (alpha * color0.g + beta * color1.g + gamma * color2.g);
                color.b = (alpha * color0.b + beta * color1.b + gamma * color2.b);
                image[x][y] = color;
            }
        }
    }
}
double calculatef(Vec3 v0, Vec3 v1, double x, double y)
{
    double x0 = v0.x;
    double x1 = v1.x;
    double y0 = v0.y;
    double y1 = v1.y;
    return (x * (y0 - y1) + y * (x1 - x0) + x0 * y1 - y0 * x1);
}

void swapVector(Vec3 &v0, Vec3 &v1)
{
    v0.x = (int)v0.x;
    v0.y = (int)v0.y;
    v0.z = (int)v0.z;
    v1.x = (int)v1.x;
    v1.y = (int)v1.y;
    v1.z = (int)v1.z;

    int temp;
    temp = (int)v0.x;
    v0.x = v1.x;
    v1.x = temp;

    temp = (int)v0.y;
    v0.y = v1.y;
    v1.y = temp;

    temp = (int)v0.z;
    v0.z = v1.z;
    v1.z = temp;

    temp = v0.colorId;
    v0.colorId = v1.colorId;
    v1.colorId = temp;
}

void midpointInterpolation(Vec3 v0, Vec3 v1)
{
    int x, y, d;

    Color c0 = colors[v0.colorId];
    Color c1 = colors[v1.colorId];
    Color color, dc;
    Color tempColor;

    // straight lines start
    if (v0.x == v1.x)
    {
        // First set color to c0
        color.r = c0.r;
        color.g = c0.g;
        color.b = c0.b;

        if (v0.y > v1.y)
        {

            dc.r = (c1.r - c0.r) / (v0.y - v1.y);
            dc.g = (c1.g - c0.g) / (v0.y - v1.y);
            dc.b = (c1.b - c0.b) / (v0.y - v1.y);

            for (int i = v0.y; i > v1.y; i--)
            {
                // round the color r,g,b values
                // then set it to image

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[(int)v0.x][i] = tempColor;

                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
        else
        {

            dc.r = (c1.r - c0.r) / (v1.y - v0.y);
            dc.g = (c1.g - c0.g) / (v1.y - v0.y);
            dc.b = (c1.b - c0.b) / (v1.y - v0.y);

            for (int i = v0.y; i < v1.y; i++)
            {
                // round the color r,g,b values
                // then set it to image

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[(int)v0.x][i] = tempColor;

                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
    }
    else if (v0.y == v1.y)
    {
        // First set color to c0
        color.r = c0.r;
        color.g = c0.g;
        color.b = c0.b;

        if (v0.x > v1.x)
        {
            dc.r = (c1.r - c0.r) / (v0.x - v1.x);
            dc.g = (c1.g - c0.g) / (v0.x - v1.x);
            dc.b = (c1.b - c0.b) / (v0.x - v1.x);

            for (int i = v0.x; i > v1.x; i--)
            {

                // round the color r,g,b values
                // then set it to image

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[i][(int)v0.y] = tempColor;

                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
        else
        {
            dc.r = (c1.r - c0.r) / (v1.x - v0.x);
            dc.g = (c1.g - c0.g) / (v1.x - v0.x);
            dc.b = (c1.b - c0.b) / (v1.x - v0.x);

            for (int i = v0.x; i < v1.x; i++)
            {
                // round the color r,g,b values
                // then set it to image

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[i][(int)v0.y] = tempColor;

                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
    }
    // straight lines end

    else if (v1.x > v0.x && v1.y < v0.y)
    {
        // 1
        // Initialize colors and color ratio changes here
        Color c0 = colors[v0.colorId];
        Color c1 = colors[v1.colorId];
        Color color, dc;

        // First set color to c0
        color.r = c0.r;
        color.g = c0.g;
        color.b = c0.b;
        if (abs(v1.y - v0.y) > abs(v1.x - v0.x))
        {

            dc.r = (c1.r - c0.r) / (v0.y - v1.y);
            dc.g = (c1.g - c0.g) / (v0.y - v1.y);
            dc.b = (c1.b - c0.b) / (v0.y - v1.y);

            x = (int)(v0.x + 0.5);
            d = 2 * (v0.x - v1.x) + (v1.y - v0.y);
            for (int y = (int)v0.y; y >= v1.y; y--)
            {
                // round the color r,g,b values
                // then set it to image
                Color tempColor;

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[x][y] = tempColor;
                if (d < 0)
                { // choose NE
                    x++;
                    d += 2 * (v0.x - v1.x + v0.y - v1.y);
                }
                else
                { // choose E
                    d += 2 * (v0.x - v1.x);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
        else
        {
            dc.r = (c1.r - c0.r) / (v1.x - v0.x);
            dc.g = (c1.g - c0.g) / (v1.x - v0.x);
            dc.b = (c1.b - c0.b) / (v1.x - v0.x);

            int y = (int)(v0.y + 0.5);
            d = 2 * (v0.y - v1.y) + (v1.x - v0.x);

            for (int x = (int)v0.x; x < v1.x; x++)
            {
                // round the color r,g,b values
                // then set it to image
                Color tempColor;

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[x][y] = tempColor;

                if (d < 0)
                { // choose NE
                    y--;
                    d += 2 * (v1.y - v0.y + v1.x - v0.x);
                }
                else
                { // choose E
                    d += 2 * (v1.y - v0.y);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
    }

    else if (v1.x < v0.x && v1.y < v0.y)
    {
        // 2

        // change v0 and v1 locations
        swapVector(v0, v1);

        // Initialize colors and color ratio changes here
        Color c0 = colors[v0.colorId];
        Color c1 = colors[v1.colorId];
        Color color, dc;

        // First set color to c0
        color.r = c0.r;
        color.g = c0.g;
        color.b = c0.b;

        if (abs(v1.y - v0.y) > abs(v1.x - v0.x))
        {
            dc.r = (c1.r - c0.r) / (v1.y - v0.y);
            dc.g = (c1.g - c0.g) / (v1.y - v0.y);
            dc.b = (c1.b - c0.b) / (v1.y - v0.y);

            x = (int)(v0.x + 0.5);
            d = 2 * (v0.x - v1.x) + (v1.y - v0.y);
            for (int y = (int)v0.y; y < v1.y - 1; y++)
            {
                // round the color r,g,b values
                // then set it to image
                Color tempColor;

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[x][y] = tempColor;

                if (d < 0)
                { // choose NE
                    x++;
                    d += 2 * (v0.x - v1.x + v1.y - v0.y);
                }
                else
                { // choose E
                    d += 2 * (v0.x - v1.x);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }

        else
        {
            dc.r = (c1.r - c0.r) / (v1.x - v0.x);
            dc.g = (c1.g - c0.g) / (v1.x - v0.x);
            dc.b = (c1.b - c0.b) / (v1.x - v0.x);

            y = (int)(v0.y + 0.5);
            d = 2 * (v0.y - v1.y) + (v1.x - v0.x);

            for (int x = (int)v0.x; x < v1.x; x++)
            {
                // round the color r,g,b values
                // then set it to image
                Color tempColor;

                tempColor.r = round(color.r);
                tempColor.g = round(color.g);
                tempColor.b = round(color.b);
                image[x][y] = tempColor;
                if (d < 0)
                { // choose NE
                    y++;
                    d += 2 * (v0.y - v1.y + v1.x - v0.x);
                }
                else
                { // choose E
                    d += 2 * (v0.y - v1.y);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
    }
    /*
    else if (v1.x < v0.x && v1.y > v0.y)
    {
        // 3

        // change v0 and v1 locations
        swapVector(v0, v1);

        // Initialize colors and color ratio changes here
        Color c0 = colors[v0.colorId];
        Color c1 = colors[v1.colorId];
        Color color, dc;

        // First set color to c0
        color.r = c0.r;
        color.g = c0.g;
        color.b = c0.b;

        if (abs(v1.y - v0.y) > abs(v1.x - v0.x))
        {
            dc.r = (c1.r - c0.r) / (v1.y - v0.y);
            dc.g = (c1.g - c0.g) / (v1.y - v0.y);
            dc.b = (c1.b - c0.b) / (v1.y - v0.y);

            x = (int)v0.x;
            d = 2 * (v0.x - v1.x) + (v1.y - v0.y);
            for (int y = (int)v0.y; y > v1.y; y--)
            {
                // round the color r,g,b values
                // then set it to image
                color.r = round(color.r);
                color.g = round(color.g);
                color.b = round(color.b);

                image[x][y] = color;
                if (d < 0)
                { // choose NE
                    x++;
                    d += 2 * (v0.x - v1.x + v0.y - v1.y);
                }
                else
                { // choose E
                    d += 2 * (v0.x - v1.x);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
        else
        {
            dc.r = (c1.r - c0.r) / (v1.x - v0.x);
            dc.g = (c1.g - c0.g) / (v1.x - v0.x);
            dc.b = (c1.b - c0.b) / (v1.x - v0.x);

            int y = (int)v0.y;
            d = 2 * (v0.y - v1.y) + (v1.x - v0.x);
            for (int x = (int)v0.x; x < v1.x; x++)
            {
                // round the color r,g,b values
                // then set it to image
                color.r = round(color.r);
                color.g = round(color.g);
                color.b = round(color.b);
                image[x][y] = color;

                if (d < 0)
                { // choose NE
                    y--;
                    d += 2 * (v1.y - v0.y + v1.x - v0.x);
                }
                else
                { // choose E
                    d += 2 * (v1.y - v0.y);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
    }

    else if (v1.x > v0.x && v1.y > v0.y)
    {
        // 4
        // Initialize colors and color ratio changes here
        Color c0 = colors[v0.colorId];
        Color c1 = colors[v1.colorId];
        Color color, dc;

        // First set color to c0
        color.r = c0.r;
        color.g = c0.g;
        color.b = c0.b;

        if (abs(v1.y - v0.y) > abs(v1.x - v0.x))
        {
            dc.r = (c1.r - c0.r) / (v1.y - v0.y);
            dc.g = (c1.g - c0.g) / (v1.y - v0.y);
            dc.b = (c1.b - c0.b) / (v1.y - v0.y);

            x = (int)v0.x;
            d = 2 * (v0.x - v1.x) + (v1.y - v0.y);
            for (int y = (int)v0.y; y < v1.y; y++)
            {
                // round the color r,g,b values
                // then set it to image
                color.r = round(color.r);
                color.g = round(color.g);
                color.b = round(color.b);

                image[x][y] = color;
                if (d < 0)
                { // choose NE
                    x++;
                    d += 2 * (v0.x - v1.x + v1.y - v0.y);
                }
                else
                { // choose E
                    d += 2 * (v0.x - v1.x);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
        else
        {
            dc.r = (c1.r - c0.r) / (v1.x - v0.x);
            dc.g = (c1.g - c0.g) / (v1.x - v0.x);
            dc.b = (c1.b - c0.b) / (v1.x - v0.x);

            y = (int)v0.y;
            d = 2 * (v0.y - v1.y) + (v1.x - v0.x);

            for (int x = (int)v0.x; x < v1.x; x++)
            {
                // round the color r,g,b values
                // then set it to image
                color.r = round(color.r);
                color.g = round(color.g);
                color.b = round(color.b);
                image[x][y] = color;

                if (d < 0)
                { // choose NE
                    y++;
                    d += 2 * (v0.y - v1.y + v1.x - v0.x);
                }
                else
                { // choose E
                    d += 2 * (v0.y - v1.y);
                }
                color.r += dc.r;
                color.g += dc.g;
                color.b += dc.b;
            }
        }
    }*/
}

int main(int argc, char **argv)
{
    int i, j;

    if (argc < 2)
    {
        std::cout << "Usage: ./rasterizer <scene file> <camera file>" << std::endl;
        return 1;
    }

    // read camera and scene files
    readSceneFile(argv[1]);
    readCameraFile(argv[2]);

    image = 0;

    for (i = 0; i < numberOfCameras; i++)
    {

        // allocate memory for image
        if (image)
        {
            for (j = 0; j < cameras[i].sizeX; j++)
            {
                delete image[j];
            }

            delete[] image;
        }

        image = new Color *[cameras[i].sizeX];

        if (image == NULL)
        {
            std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
            exit(1);
        }

        for (j = 0; j < cameras[i].sizeX; j++)
        {
            image[j] = new Color[cameras[i].sizeY];
            if (image[j] == NULL)
            {
                std::cout << "ERROR: Cannot allocate memory for image." << std::endl;
                exit(1);
            }
        }

        // initialize image with basic values
        initializeImage(cameras[i]);

        // do forward rendering pipeline operations
        forwardRenderingPipeline(cameras[i]);

        // generate PPM file
        writeImageToPPMFile(cameras[i]);

        // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
        // Notice that os_type is not given as 1 (Ubuntu) or 2 (Windows), below call doesn't do conversion.
        // Change os_type to 1 or 2, after being sure that you have ImageMagick installed.
        convertPPMToPNG(cameras[i].outputFileName, 99);
    }

    return 0;
}