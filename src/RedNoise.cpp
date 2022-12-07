#include "../libs/sdw/CanvasTriangle.h"
#include "../libs/sdw/DrawingWindow.h"
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <algorithm>
#include <Utils.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <glm/glm.hpp>
#include "glm/ext.hpp"
#include <CanvasPoint.h>
#include <Colour.h>
#include <cmath>

//h = 320
#define WIDTH 640
#define HEIGHT 480

glm::vec3 cam = glm::vec3(0.0, 0.0, 8.0);
float focalLength = 2.0;
glm::vec3 lightSource = glm::vec3(0.0, 2.0, 0.0);
glm::mat3 camOrientation = glm::mat3(-1, 0, 0,
                                     0, 1, 0,
                                     0, 0, 1);
int orbit;
int mode;
glm::vec3 camTest = glm::vec3(0.0, 0.0, 4.0);
//float depth[HEIGHT][WIDTH];


std::vector<std::vector<float>> depth(HEIGHT, std::vector<float>(WIDTH));
//std::vector<std::vector<float>> colourBuffer(HEIGHT, std::vector<float>(WIDTH));

void printVector(std::vector<float> input, std::string name){
    std::cout << name << ":";
    for (int i = 0; i < input.size(); ++i) {
        std::cout << input.at(i) << " ";
    }
    std::cout << std::endl;
}

std::vector<float> interpolateSingleFloats(float from, float to, float numberOfValues){
    if (numberOfValues < 0 ){
        numberOfValues *= -1;
    }
    float space;
    std::vector<float> output;
    if ((int)numberOfValues == 0){
        output.push_back(from);
        output.push_back(to);
    } else {
        space = (to - from)/(float)abs(numberOfValues);
        for (int i = 0; i < numberOfValues+1; ++i) {
            float v = from + (space * (float)i);
            output.push_back(v);
        }
    }


    return output;
}
//task 4
std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues){
    //first column
    std::vector<float> fCol = interpolateSingleFloats(from[0], to[0], numberOfValues);
    //second column
    std::vector<float> sCol = interpolateSingleFloats(from[1], to[1], numberOfValues);
    //third column
    std::vector<float> tCol = interpolateSingleFloats(from[2], to[2], numberOfValues);

    std::vector<glm::vec3> result;
    for (int i = 0; i < numberOfValues; ++i) {

        std::vector<float> copy;
        copy.push_back(fCol[i]);
        copy.push_back(sCol[i]);
        copy.push_back(tCol[i]);

        glm::vec3 toAppend(copy[0],copy[1],copy[2]);
        result.push_back(toAppend);

        copy.clear();
    }
    return result;
}

//// CAMERA AND LIGHT FUNCTIONS

void getLight(glm::vec3 input, glm::vec3 offset){
    lightSource = input + offset;
}

void rotateLight(glm::vec3 input, glm::mat3 offset){
    lightSource = offset * input;
}

void getCamera(glm::vec3 input, glm::vec3 offset){
    cam = input + offset;
}

void getFocal(float input, float offset){
    focalLength = input + offset;
}

void rotateCamera(glm::vec3 input, glm::mat3 offset){
    cam = offset * input;
}

void resetOrientation(){
    camOrientation = glm::mat3(1, 0, 0,
                               0, -1, 0,
                               0, 0, 1);
}

void lookAt(glm::vec3 target){
    resetOrientation();
    glm::vec3 forward, up, right, to;
    glm::vec3 vertical = glm::vec3(0, -1, 0);
    to = target - cam;
    forward = glm::normalize(cam - to);
    right = glm::normalize(glm::cross(vertical, forward));
    up = glm::normalize(glm::cross(forward, right));


    glm::mat3 mat = glm::mat3(right.x, up.x, forward.x,
                              right.y, up.y, forward.y,
                              right.z, up.z, forward.z);

    camOrientation = camOrientation * mat;
}

void setDepthToZero(){
    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
            depth[y][x] = 0;
        }
    }
}

//// POINTS CALCULATION FUNCTIONS

std::vector<float> interpolatePoints(CanvasPoint from, CanvasPoint to){

    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff/numberOfSteps;
    float yStepSize = yDiff/numberOfSteps;

    std::vector<float> output;
    output.push_back(from.x);
    output.push_back(from.y);
    output.push_back(xStepSize);
    output.push_back(yStepSize);
    output.push_back(numberOfSteps);
    return output;
}
//task 2
void drawLine(DrawingWindow &window, CanvasPoint start, CanvasPoint end, float leftZ, float rightZ, uint32_t colour){

    std::vector<float> stepSizes = interpolatePoints(start, end);
    std::vector<float> zDepth = interpolateSingleFloats(leftZ, rightZ, stepSizes[4]);
    std::vector<float> xCoord = interpolateSingleFloats(start.x, end.x, stepSizes[4]);

//    printVector(zDepth, "zDepth");
//    std::cout << (int)stepSizes[4] << std::endl;

    if (leftZ == 0 && rightZ == 0){
        for (int i = 0; i < stepSizes[4]; ++i) {
            float x = stepSizes[0] + (stepSizes[2]*(float)i);
            float y = stepSizes[1] + (stepSizes[3]*(float)i);
            if (x > WIDTH || x < 0 || y > HEIGHT || y < 0){
                continue;
            } else {
                window.setPixelColour((int) round(x), (int) round(y), colour);
            }
        }
    } else {
        for (float i = 0; i < stepSizes[4]; ++i) {
            float x = stepSizes[0] + (stepSizes[2]*(float)i);
            float y = start.y;
            if (x > WIDTH || x < 0 || x < start.x || x > end.x || y > HEIGHT || y < 0){
                continue;
            }
//            if (x >= start.x && x <= end.x && x >= 0 && x <= WIDTH){
            if (depth[(int)round(y)][(int)round(x)] == 0){
                if (round(x) >= start.x && round(x) <= end.x){
                    if (x > WIDTH || x < 0 || x < start.x || x > end.x || y > HEIGHT || y < 0){
                        continue;
                    } else {
                        window.setPixelColour((int) round(x), (int) round(y), colour);
                        depth[(int)round(y)][(int)round(x)] = zDepth[(int)round(i)];
                    }
                }
            } else {
                if (zDepth[(int)round(i)] >= depth[(int)round(y)][(int)round(x)]){
                    if (round(x) >= start.x && round(x) <= end.x) {
                        if (x > WIDTH || x < 0 || x < start.x || x > end.x || y > HEIGHT || y < 0){
                            continue;
                        } else {
                            window.setPixelColour((int) round(x), (int) round(y), colour);
                            depth[(int)round(y)][(int)round(x)] = zDepth[(int)round(i)];
                        }
                    }
                }else continue;
            }
//            }
        }
    }

}

void drawTexture(DrawingWindow &window, CanvasPoint start, CanvasPoint end, std::vector<std::vector<uint32_t>> matrix){
//    window.clearPixels();

    //stepSizes(originX, originY, xStepSize, yStepSize, numberOfSteps)
    std::vector<float> stepSizes = interpolatePoints(start, end);

//    uint32_t colour = (255 << 24) + (rand()%255 << 16) + (rand()%255 << 8) + rand()%255 ;
//    uint32_t colour = uint32_t(start.depth) ;

    for (float i = 0; i < stepSizes[4]; ++i) {
        float x = stepSizes[0] + (stepSizes[2]*i);
        float y = stepSizes[1] + (stepSizes[3]*i);
        window.setPixelColour(round(x), round(y), matrix[round(y)][round(x)]);
    }

}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, uint32_t colour){
//    window.clearPixels();

    std::vector<float> side1 = interpolateSingleFloats(triangle.v1().depth, triangle.v0().depth, fabs(triangle.v1().x - triangle.v0().x));
    std::vector<float> side2 = interpolateSingleFloats(triangle.v2().depth, triangle.v1().depth, fabs(triangle.v2().x - triangle.v1().x));
    std::vector<float> side3 = interpolateSingleFloats(triangle.v0().depth, triangle.v2().depth, fabs(triangle.v0().x - triangle.v2().x));


    drawLine(window, triangle.v0(), triangle.v1(), 0, 0, colour);
    drawLine(window, triangle.v1(), triangle.v2(), 0, 0, colour);
    drawLine(window, triangle.v2(), triangle.v0(), 0, 0, colour);
}

float findPointX(std::vector<CanvasPoint> line, float pointY){
    float pointX;
    float gradient;
    float yDiff = line[1].y - line[0].y;
    float xDiff = line[1].x - line[0].x;
    if (xDiff == 0){
        gradient = 0;
        pointX = line[0].x;
    } else {
        gradient = yDiff/xDiff;
        float c = line[0].y - (gradient*line[0].x);
        pointX = (pointY - c)/gradient;
    }

    return pointX;
}

//// RASTERIZER
void fillTriangle(DrawingWindow &window, CanvasTriangle triangle, uint32_t fillColour){

    CanvasPoint mid;
    std::vector<CanvasPoint> line;  //adjacent line
//    drawTriangle(window, triangle, fillColour);

    float maxY = std::max(triangle.v0().y, std::max(triangle.v1().y, triangle.v2().y));
    float minY = std::min(triangle.v0().y, std::min(triangle.v1().y, triangle.v2().y));
    float maxX = std::max(triangle.v0().x, std::max(triangle.v1().x, triangle.v2().x));
    float minX = std::min(triangle.v0().x, std::min(triangle.v1().x, triangle.v2().x));

    // find middle point
    // first case if triangle has a flat base
    if (triangle.v0().y == triangle.v1().y || triangle.v0().y == triangle.v2().y || triangle.v1().y == triangle.v2().y){
//        std::cout << "flat based" << std::endl;
        CanvasPoint leftPoint;
        CanvasPoint rightPoint;
        CanvasPoint otherPoint;

        std::vector<float> leftSide;
        std::vector<float> rightSide;
//        std::vector<float> zDepths;

        if (triangle.v0().y == triangle.v1().y) {
            float left = std::min(triangle.v0().x, triangle.v1().x);

            otherPoint = triangle.v2();
            if (left == triangle.v0().x){
                leftPoint = triangle.v0();
                rightPoint = triangle.v1();
            } else if (left == triangle.v1().x){
                leftPoint = triangle.v1();
                rightPoint = triangle.v0();
            }
        } else if (triangle.v0().y == triangle.v2().y) {
            float left = std::min(triangle.v0().x, triangle.v2().x);

            otherPoint = triangle.v1();
            if (left == triangle.v0().x){
                leftPoint = triangle.v0();
                rightPoint = triangle.v2();
            } else if (left == triangle.v2().x){
                leftPoint = triangle.v2();
                rightPoint = triangle.v0();
            }
        } else if (triangle.v1().y == triangle.v2().y) {
            float left = std::min(triangle.v1().x, triangle.v2().x);

            otherPoint = triangle.v0();
            if (left == triangle.v1().x){
                leftPoint = triangle.v1();
                rightPoint = triangle.v2();
            } else if (left == triangle.v2().x){
                leftPoint = triangle.v2();
                rightPoint = triangle.v1();
            }
        }

        // interpolate zDepth of left and right border
        if (otherPoint.y < leftPoint.y){
            leftSide = interpolateSingleFloats(otherPoint.depth, leftPoint.depth, (otherPoint.y - leftPoint.y));
            rightSide = interpolateSingleFloats(otherPoint.depth, rightPoint.depth, (otherPoint.y - rightPoint.y));
        } else if (otherPoint.y > leftPoint.y){
            leftSide = interpolateSingleFloats(leftPoint.depth, otherPoint.depth, (otherPoint.y - leftPoint.y));
            rightSide = interpolateSingleFloats(rightPoint.depth, otherPoint.depth, (otherPoint.y - rightPoint.y));
        }

//        printVector(leftSide, "leftSide");
//        printVector(rightSide, "rightSide");

        int index = 0;
        for (float y = minY; y <= maxY; ++y) {
            if (y == otherPoint.y){
                if (otherPoint.x > WIDTH || otherPoint.x < 0 || otherPoint.y > HEIGHT || otherPoint.y < 0){
                    continue;
                } else {
                    if (depth[(int)otherPoint.y][(int)otherPoint.x] == 0){
                        window.setPixelColour((int)round(otherPoint.x), (int)round(otherPoint.y), fillColour);
                        depth[(int)round(otherPoint.y)][(int)round(otherPoint.x)] = otherPoint.depth;
                    } else if (otherPoint.depth > depth[(int)otherPoint.y][(int)otherPoint.x] && depth[(int)otherPoint.y][(int)otherPoint.x] != 0){
                        window.setPixelColour((int)round(otherPoint.x), (int)round(otherPoint.y), fillColour);
                        depth[(int)round(otherPoint.y)][(int)round(otherPoint.x)] = otherPoint.depth;
                    }
                }
            } else {
                drawLine(window,
                         CanvasPoint{findPointX(std::vector<CanvasPoint>{otherPoint, leftPoint}, y), y},
                         CanvasPoint{findPointX(std::vector<CanvasPoint>{otherPoint, rightPoint}, y), y},
                         leftSide[index], rightSide[index], fillColour);
            }
            index++;
        }
        index = 0;
    } else {
        //if triangle doesn't have a flat base
//        std::cout << "not flat" << std::endl;
        if ((triangle.v1().y < triangle.v0().y && triangle.v0().y < triangle.v2().y) ||
            (triangle.v2().y < triangle.v0().y && triangle.v0().y < triangle.v1().y)) {
            mid = triangle.v0();
            line.push_back(triangle.v1());
            line.push_back(triangle.v2());
        } else if ((triangle.v0().y < triangle.v1().y && triangle.v1().y < triangle.v2().y) ||
                   (triangle.v2().y < triangle.v1().y && triangle.v1().y < triangle.v0().y)) {
            mid = triangle.v1();
            line.push_back(triangle.v0());
            line.push_back(triangle.v2());
        } else if ((triangle.v0().y < triangle.v2().y && triangle.v2().y < triangle.v1().y) ||
                   (triangle.v1().y < triangle.v2().y && triangle.v2().y < triangle.v0().y)) {
            mid = triangle.v2();
            line.push_back(triangle.v0());
            line.push_back(triangle.v1());
        }

        //to find uppermost, lowermost, leftmost and rightmost point
        CanvasPoint top, bottom;
        float topX, bottomX;

        if (minY == triangle.v0().y){
            top = triangle.v0();
        } else if (minY == triangle.v1().y){
            top = triangle.v1();
        } else if (minY == triangle.v2().y){
            top = triangle.v2();
        }

        if (maxY == triangle.v0().y){
            bottom = triangle.v0();
        } else if (maxY == triangle.v1().y){
            bottom = triangle.v1();
        } else if (maxY == triangle.v2().y){
            bottom = triangle.v2();
        }

        std::vector<float> side1 = interpolateSingleFloats(top.depth, mid.depth, (mid.y - top.y));
        std::vector<float> side2 = interpolateSingleFloats(mid.depth, bottom.depth, (bottom.y - mid.y));
        std::vector<float> adjacentSide = interpolateSingleFloats(top.depth, bottom.depth, (bottom.y - top.y));
//        std::cout << "values:  " << (int)(mid.y - top.y) << " ";
//        printVector(side1, "side1");
//        std::cout << "values:  " << (int)(bottom.y - mid.y) << " ";
//        printVector(side2, "side2");
//        std::cout << "values:  " << (int)(bottom.y - top.y) << " ";
//        printVector(adjacentSide, "adjacentSide"); // some interpolation didn't have contents NEED TO FIX

        //draw triangle border
        float xIntercept = findPointX(line, mid.y);
//        drawTriangle(window, triangle, fillColour);
        CanvasPoint intercept = CanvasPoint{xIntercept, mid.y, adjacentSide[mid.y-top.y]};

        //fill in upper half
        int index = 0;
//        std::cout << "filling upper half" << std::endl;
        for (float y = minY; y < mid.y+1; ++y) {
            if (y == minY) {
                if (top.x > WIDTH || top.x < 0 || top.y > HEIGHT || top.y < 0){
                    continue;
                } else {
                    if (depth[(int)round(top.y)][(int)round(top.x)] == 0){
                        window.setPixelColour((int)round(top.x),(int)round(top.y), fillColour);
                        depth[(int)round(top.y)][(int)round(top.x)] = top.depth;
                    } else {
                        if (top.depth > depth[(int)round(top.y)][(int)round(top.x)]){
                            window.setPixelColour((int)round(top.x),(int)round(top.y), fillColour);
                            depth[(int)round(top.y)][(int)round(top.x)] = top.depth;
                        }
                    }
                }
            } else {
                if (mid.x < xIntercept){ //if midpoint on left side
                    drawLine(window,
                             CanvasPoint{findPointX(std::vector<CanvasPoint>{top, mid}, y), y},
                             CanvasPoint{findPointX(line, y), y},
                             side1[index], adjacentSide[index], fillColour); //segfault when interpolating
                } else if (mid.x > xIntercept){ //if midpoint on the right
                    drawLine(window,
                             CanvasPoint{findPointX(line, y), y},
                             CanvasPoint{findPointX(std::vector<CanvasPoint>{top, mid}, y), y},
                             adjacentSide[index], side1[index], fillColour);
                }
            }
            index++;
        }

        int offset = index + 1;
        index = 0;

        //fill in lower half
//        std::cout << "filling lower half" << std::endl;

        for (float y = mid.y+1; y <= maxY; ++y) {
            if (y == maxY) {
                if (bottom.x > WIDTH || bottom.x < 0 || bottom.y > HEIGHT || bottom.y < 0){
                    continue;
                } else {
                    if (depth[(int)round(bottom.y)][(int)round(bottom.x)] == 0){
                        window.setPixelColour((int)round(bottom.x),(int)round(bottom.y), fillColour);
                        depth[(int)round(bottom.y)][(int)round(bottom.x)] = bottom.depth;
                    } else {
                        if (bottom.depth >= depth[(int)round(bottom.y)][(int)round(bottom.x)]){
                            window.setPixelColour((int)round(bottom.x),(int)round(bottom.y), fillColour);
                            depth[(int)round(bottom.y)][(int)round(bottom.x)] = bottom.depth;
                        }
                    }
                }
            } else {
                if (mid.x < xIntercept){
                    drawLine(window,
                         CanvasPoint{findPointX(std::vector<CanvasPoint>{mid, bottom}, y), y},
                         CanvasPoint{findPointX(line, y), y},
                         side2[index], adjacentSide[offset+index], fillColour);
                } else if (mid.x > xIntercept){
                    drawLine(window,
                         CanvasPoint{findPointX(line, y), y},
                         CanvasPoint{findPointX(std::vector<CanvasPoint>{mid, bottom}, y), y},
                         adjacentSide[offset+index], side2[index], fillColour);
                }
                index++;
            }

        }
    }
}

//// RAY TRACING
RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPosition, glm::vec3 rayDirection, std::vector<ModelTriangle> models){
    RayTriangleIntersection output;
    float t_buffer = INFINITY;
    for (int i = 0; i < models.size(); ++i) {
        ModelTriangle triangle = models[i];

        glm::vec3 p0 = triangle.vertices[0];
        glm::vec3 p1 = triangle.vertices[1];
        glm::vec3 p2 = triangle.vertices[2];

        glm::vec3 e0 = p1 - p0;
        glm::vec3 e1 = p2 - p0;
        glm::vec3 SPVector = cam - p0; //cam -p0
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        glm::vec3 intersection_p = triangle.vertices[0] + (e0 * possibleSolution[1]) + (e1 * possibleSolution[2]);

        float t = possibleSolution[0];
        float u = possibleSolution[1];
        float v = possibleSolution[2];

        if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && t >= 0 && t < t_buffer){
            t_buffer = t;
            output = RayTriangleIntersection(intersection_p, t, triangle, 1);
        }
    }
    return output;
}

bool checkShadow(RayTriangleIntersection intersection, glm::vec3 startPoint, glm::vec3 rayDirection, std::vector<ModelTriangle> models){
    bool output = false;
    glm::vec3 r = lightSource - startPoint;
    float max_length = sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));

    for (int i = 0; i < models.size(); ++i) {

        ModelTriangle triangle = models[i];

        glm::vec3 p0 = triangle.vertices[0];
        glm::vec3 p1 = triangle.vertices[1];
        glm::vec3 p2 = triangle.vertices[2];

        glm::vec3 i0 = intersection.intersectedTriangle.vertices[0];
        glm::vec3 i1 = intersection.intersectedTriangle.vertices[1];
        glm::vec3 i2 = intersection.intersectedTriangle.vertices[2];

        if (p0 == i0 && p1 == i1 && p2 ==i2){
//            std::cout << "continue" << std::endl;
            continue;
        }

        glm::vec3 e0 = p1 - p0;
        glm::vec3 e1 = p2 - p0;
        glm::vec3 SPVector = startPoint - p0; //cam -p0
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        glm::vec3 intersection_p = triangle.vertices[0] + (e0 * possibleSolution[1]) + (e1 * possibleSolution[2]);

        float t = possibleSolution[0];
        float u = possibleSolution[1];
        float v = possibleSolution[2];

//        std::cout << "current t: " << t << std::endl;

        if ((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && t > 0){
//            output = true;
            if (t < max_length){
                return true;
            }
//            break;
        } else {
            output = false;
        }
    }
    return output;
}

//// CALCULATION FOR RENDERING

ModelTriangle calcNormal(ModelTriangle model){
    glm::vec3 vec1 = model.vertices[1] - model.vertices[0];
    glm::vec3 vec2 = model.vertices[2] - model.vertices[0];

    glm::vec3 norm = glm::normalize(glm::cross(vec1, vec2));
//    float ang = glm::dot(norm, rayDir);
//    std::cout << "normal: " << glm::to_string(norm) << std::endl;

    ModelTriangle output = ModelTriangle(model.vertices[0], model.vertices[1], model.vertices[2], model.colour, norm);

    return output;
}

glm::vec3 computeSurfaceToLight(glm::vec3 surface, glm::vec3 light){
    glm::vec3 direction = light - surface;
    return direction;
}

float surfaceToLightLength(glm::vec3 ray){
    float length = sqrt(pow(ray[0], 2) + pow(ray[1], 2) + pow(ray[2], 2));
    return length;
}

glm::vec3 computeRay(int x, int y){

    float inverseFocal = -1/(focalLength * 80);
    glm::mat3 camOrient = glm::inverse(camOrientation);
    float x_c = (x - (WIDTH/2)) * inverseFocal;
    float y_c = (y - (HEIGHT/2)) * inverseFocal;
    float z_c = -1;
    glm::vec3 rayDirection = glm::normalize(camOrient * glm::vec3(x_c, y_c, z_c));
//    std::cout << "ray_dir: " << glm::to_string(rayDirection) << std::endl;
    return rayDirection;
}

// INTERSECTION POINT ON SCREEN
CanvasPoint getCanvasIntersectionPoint(glm::vec3 cameraPosition, glm::vec3 vertexPosition, float focalLength){
    CanvasPoint point;
    float xOutput;
    float yOutput;
//    glm::mat4 view = lookAt(cameraPosition, glm::vec3(0, 0, 0));
//    lookAt(glm::vec3(0, 0, 0));
    glm::vec3 vertex = vertexPosition - cameraPosition;
    glm::vec3 vertexConverted = vertex * camOrientation;

//    std::cout << "vertex: " << glm::to_string(vertex) << std::endl;
//    std::cout << "converted: " << glm::to_string(vertexConverted) << std::endl;

    xOutput = (focalLength * (vertexConverted[0] / vertexConverted[2]) * 80 ) + (float(WIDTH) / 2);
    yOutput = (focalLength * (vertexConverted[1] / vertexConverted[2]) * 80 ) + (float(HEIGHT) / 2);
    float pointDepth  = -1 / vertexConverted[2];


    point = CanvasPoint(xOutput, yOutput, pointDepth); //depth already in camera space

    return point;
}

uint32_t getColourShaded (float length, float dp, std::vector<Colour> colourList, const ModelTriangle& model){
    uint32_t rgb;
//    int intensity = ;
    for (int i = 0; i < colourList.size(); ++i) {
        if (model.colour.name == colourList[i].name) {
            int red = int((colourList[i].red / (length/3)) * dp);
            int blue = int((colourList[i].blue / (length/3)) * dp);
            int green = int((colourList[i].green / (length/3)) * dp);
            if (red > 255){
                red = 255;
            }
            if (blue > 255){
                blue = 255;
            }
            if (green > 255){
                green = 255;
            }
            rgb = (255 << 24) + (red << 16) + (green << 8) + blue;
        }
    }
    return rgb;
}

uint32_t getColourShadow (float length, std::vector<Colour> colourList, const ModelTriangle& model){
    uint32_t rgb;
//    int intensity = ;
    for (int i = 0; i < colourList.size(); ++i) {
        if (model.colour.name == colourList[i].name) {
            int red = int(colourList[i].red / length);
            int blue = int(colourList[i].blue / length);
            int green = int(colourList[i].green / length);
            int alpha = 0;
            if (red > 255){
                red = 255;
            }
            if (blue > 255){
                blue = 255;
            }
            if (green > 255){
                green = 255;
            }
            rgb = (alpha << 24) + (red << 16) + (green << 8) + blue;
        }
    }
    return rgb;
}

uint32_t getColour (std::vector<Colour> colourList, const ModelTriangle& model){
    uint32_t rgb;
    for (int i = 0; i < colourList.size(); ++i) {
        if (model.colour.name == colourList[i].name) {
//                rgb = (0 << 24) + (colourList[i].red << 16) + (colourList[i].green << 8) + colourList[i].blue;
//                std::cout << "shadow: " << rgb << std::endl;
                rgb = (255 << 24) + (colourList[i].red << 16) + (colourList[i].green << 8) + colourList[i].blue;
//                std::cout << "not shadow: " << rgb << std::endl;
        }
    }
    return rgb;
}

void drawRaster(DrawingWindow &window, const std::string& objfile, const std::string& mtlfile, glm::vec3 cameraPosition){

    //READ .mtl FILE
    std::ifstream mtlFile (mtlfile);
    std::string mtlLine;
    std::string newmtl = "newmtl";
    std::string kd = "Kd";
    std::vector<std::string> colourName;
    std::vector<std::vector<float>> rgbColours;
    std::vector<Colour> colourList;
    uint32_t colourConverted;

    while(!mtlFile.eof()){
        getline(mtlFile, mtlLine);
        int search1 = mtlLine.find(newmtl);
        int search2 = mtlLine.find(kd);
        if (search1 != std::string::npos){
            std::string tempColour;
            mtlLine.erase(search1, newmtl.length()+1);
//            std::cout << mtlLine << std::endl;
            std::istringstream iss(mtlLine);
            iss >> tempColour;
            colourName.push_back(tempColour);
        } else if (search2 != std::string::npos){
            std::vector<float> rgb(3);
            mtlLine.erase(search2, kd.length()+1);
//            std::cout << mtlLine << std::endl;
            std::istringstream iss(mtlLine);
            iss >> rgb[0] >> rgb[1] >> rgb[2];
            rgbColours.push_back(rgb);
        }
    }

    //Organize colours
    for (int i = 0; i < colourName.size(); ++i) {
        Colour tempColour;
        tempColour = Colour(colourName[i],
                            int(rgbColours[i][0]*255),
                            int(rgbColours[i][1]*255),
                            int(rgbColours[i][2]*255));
        colourList.push_back(tempColour);

//        std::cout << tempColour << std::endl;
    }

    //READ .obj FILE
    // process 'v' coordinates
    ModelTriangle model;
    std::vector<ModelTriangle> coordinates;
    std::string v, valuesX, valuesY, valuesZ;
    std::vector<glm::vec3> points;
//    std::vector<std::vector<CanvasPoint>> toFill; // coordinates with depth data
    // process 'f' (facets)
    std::string f, f1, f2, f3;
    std::vector<std::vector<int>> facets;

    std::string line;
    std::ifstream objFile (objfile);
    std::string usemtl = "usemtl";
    Colour colour;

    while (!objFile.eof()){
        getline(objFile, line);
        int searchColour = line.find(usemtl);
        if (line[0] == 'v'){
            std::istringstream iss(line);
//            std::cout << line << std::endl;
            iss >> v >> valuesX >> valuesY >> valuesZ;
            glm::vec3 point(std::stof(valuesX) * 0.35, std::stof(valuesY) * 0.35, std::stof(valuesZ) * 0.35);
            points.push_back(point);
            v.clear(); valuesX.clear(); valuesY.clear(); valuesZ.clear();
        } else if (line[0] == 'f'){
            // remove "/" in each line to extract the numbers
            line.erase(std::remove(line.begin(), line.end(), '/'), line.end());
//            std::cout << line << std::endl;
            std::istringstream iss(line);
            iss >> f >> f1 >> f2 >> f3;
            // build model based on facets
            model = ModelTriangle(points[std::stoi(f1)-1],
                                  points[std::stoi(f2)-1],
                                  points[std::stoi(f3)-1],
                                  colour );
            f.clear(); f1.clear(); f2.clear(); f3.clear();
//            std::cout << model << std::endl;

            //// render directly
            std::vector<CanvasPoint> temp;
            for (int j = 0; j < 3; ++j) {
                CanvasPoint screen = getCanvasIntersectionPoint(cameraPosition,
                                                                model.vertices[j],
                                                                focalLength);
                temp.push_back(screen);
            }
//            std::cout << model.colour.name << std::endl;
            fillTriangle(window, CanvasTriangle(temp[0], temp[1], temp[2]), getColour(colourList, model));
            temp.clear();
            coordinates.push_back(model);

        } else if(searchColour != std::string::npos){
            std::string tempColour;
            line.erase(searchColour, usemtl.length()+1);
            std::istringstream iss(line);
//            std::cout << line << std::endl;
            iss >> tempColour;
            //set colour everytime encounter with 'usemtl'
            for (int i = 0; i < colourList.size(); ++i) {
                if (tempColour == colourList[i].name){
                    colour = colourList[i];
                }
            }
        }
    }

}

void drawWireframe(DrawingWindow &window, const std::string& objfile, const std::string& mtlfile, glm::vec3 cameraPosition){

    //READ .mtl FILE
    std::ifstream mtlFile (mtlfile);
    std::string mtlLine;
    std::string newmtl = "newmtl";
    std::string kd = "Kd";
    std::vector<std::string> colourName;
    std::vector<std::vector<float>> rgbColours;
    std::vector<Colour> colourList;
    uint32_t colourConverted;

    while(!mtlFile.eof()){
        getline(mtlFile, mtlLine);
        int search1 = mtlLine.find(newmtl);
        int search2 = mtlLine.find(kd);
        if (search1 != std::string::npos){
            std::string tempColour;
            mtlLine.erase(search1, newmtl.length()+1);
//            std::cout << mtlLine << std::endl;
            std::istringstream iss(mtlLine);
            iss >> tempColour;
            colourName.push_back(tempColour);
        } else if (search2 != std::string::npos){
            std::vector<float> rgb(3);
            mtlLine.erase(search2, kd.length()+1);
//            std::cout << mtlLine << std::endl;
            std::istringstream iss(mtlLine);
            iss >> rgb[0] >> rgb[1] >> rgb[2];
            rgbColours.push_back(rgb);
        }
    }

    //Organize colours
    for (int i = 0; i < colourName.size(); ++i) {
        Colour tempColour;
        tempColour = Colour(colourName[i],
                            int(rgbColours[i][0]*255),
                            int(rgbColours[i][1]*255),
                            int(rgbColours[i][2]*255));
        colourList.push_back(tempColour);

//        std::cout << tempColour << std::endl;
    }

    //READ .obj FILE
    // process 'v' coordinates
    ModelTriangle model;
    std::vector<ModelTriangle> coordinates;
    std::string v, valuesX, valuesY, valuesZ;
    std::vector<glm::vec3> points;
//    std::vector<std::vector<CanvasPoint>> toFill; // coordinates with depth data
    // process 'f' (facets)
    std::string f, f1, f2, f3;
    std::vector<std::vector<int>> facets;

    std::string line;
    std::ifstream objFile (objfile);
    std::string usemtl = "usemtl";
    Colour colour;

    while (!objFile.eof()){
        getline(objFile, line);
        int searchColour = line.find(usemtl);
        if (line[0] == 'v'){
            std::istringstream iss(line);
//            std::cout << line << std::endl;
            iss >> v >> valuesX >> valuesY >> valuesZ;
            glm::vec3 point(std::stof(valuesX) * 0.35, std::stof(valuesY) * 0.35, std::stof(valuesZ) * 0.35);
            points.push_back(point);
            v.clear(); valuesX.clear(); valuesY.clear(); valuesZ.clear();
        } else if (line[0] == 'f'){
            // remove "/" in each line to extract the numbers
            line.erase(std::remove(line.begin(), line.end(), '/'), line.end());
//            std::cout << line << std::endl;
            std::istringstream iss(line);
            iss >> f >> f1 >> f2 >> f3;
            // build model based on facets
            model = ModelTriangle(points[std::stoi(f1)-1],
                                  points[std::stoi(f2)-1],
                                  points[std::stoi(f3)-1],
                                  colour );
            f.clear(); f1.clear(); f2.clear(); f3.clear();
//            std::cout << model << std::endl;

            //// render directly
            std::vector<CanvasPoint> temp;
            for (int j = 0; j < 3; ++j) {
                CanvasPoint screen = getCanvasIntersectionPoint(cameraPosition,
                                                                model.vertices[j],
                                                                focalLength);
                temp.push_back(screen);
            }
//            std::cout << model.colour.name << std::endl;
            drawTriangle(window, CanvasTriangle(temp[0], temp[1], temp[2]), getColour(colourList, model));
            temp.clear();
            coordinates.push_back(model);

        } else if(searchColour != std::string::npos){
            std::string tempColour;
            line.erase(searchColour, usemtl.length()+1);
            std::istringstream iss(line);
//            std::cout << line << std::endl;
            iss >> tempColour;
            //set colour everytime encounter with 'usemtl'
            for (int i = 0; i < colourList.size(); ++i) {
                if (tempColour == colourList[i].name){
                    colour = colourList[i];
                }
            }
        }
    }

}

void drawRayTrace(DrawingWindow &window, const std::string& objfile, const std::string& mtlfile, glm::vec3 cameraPosition){

    //READ .mtl FILE
    std::ifstream mtlFile (mtlfile);
    std::string mtlLine;
    std::string newmtl = "newmtl";
    std::string kd = "Kd";
    std::vector<std::string> colourName;
    std::vector<std::vector<float>> rgbColours;
    std::vector<Colour> colourList;
    uint32_t colourConverted;

    while(!mtlFile.eof()){
        getline(mtlFile, mtlLine);
        int search1 = mtlLine.find(newmtl);
        int search2 = mtlLine.find(kd);
        if (search1 != std::string::npos){
            std::string tempColour;
            mtlLine.erase(search1, newmtl.length()+1);
//            std::cout << mtlLine << std::endl;
            std::istringstream iss(mtlLine);
            iss >> tempColour;
            colourName.push_back(tempColour);
        } else if (search2 != std::string::npos){
            std::vector<float> rgb(3);
            mtlLine.erase(search2, kd.length()+1);
//            std::cout << mtlLine << std::endl;
            std::istringstream iss(mtlLine);
            iss >> rgb[0] >> rgb[1] >> rgb[2];
            rgbColours.push_back(rgb);
        }
    }

    //Organize colours
    for (int i = 0; i < colourName.size(); ++i) {
        Colour tempColour;
        tempColour = Colour(colourName[i],
                            int(rgbColours[i][0]*255),
                            int(rgbColours[i][1]*255),
                            int(rgbColours[i][2]*255));
        colourList.push_back(tempColour);

//        std::cout << tempColour << std::endl;
    }

    //READ .obj FILE
    // process 'v' coordinates
    std::vector<ModelTriangle> coordinates;
    std::string v, valuesX, valuesY, valuesZ;
    std::vector<glm::vec3> points;
//    std::vector<std::vector<CanvasPoint>> toFill; // coordinates with depth data
    // process 'f' (facets)
    std::string f, f1, f2, f3;
    std::vector<std::vector<int>> facets;

    std::string line;
    std::ifstream objFile (objfile);
    std::string usemtl = "usemtl";
    Colour colour;

    while (!objFile.eof()){
        getline(objFile, line);
        int searchColour = line.find(usemtl);
        if (line[0] == 'v'){
            std::istringstream iss(line);
//            std::cout << line << std::endl;
            iss >> v >> valuesX >> valuesY >> valuesZ;
            glm::vec3 point(std::stof(valuesX), std::stof(valuesY), std::stof(valuesZ));
            points.push_back(point);
            v.clear(); valuesX.clear(); valuesY.clear(); valuesZ.clear();
        } else if (line[0] == 'f'){
            // remove "/" in each line to extract the numbers
            line.erase(std::remove(line.begin(), line.end(), '/'), line.end());
//            std::cout << line << std::endl;
            std::istringstream iss(line);
            iss >> f >> f1 >> f2 >> f3;
            // build model based on facets
            ModelTriangle model_n = ModelTriangle(points[std::stoi(f1)-1],
                                                points[std::stoi(f2)-1],
                                                points[std::stoi(f3)-1],
                                                colour);
            ModelTriangle model = calcNormal(model_n);
//            std::cout << "normal: " << glm::to_string(model.normal) << std::endl;
            f.clear(); f1.clear(); f2.clear(); f3.clear();

            coordinates.push_back(model);

        } else if(searchColour != std::string::npos){
            std::string tempColour;
            line.erase(searchColour, usemtl.length()+1);
            std::istringstream iss(line);
//            std::cout << line << std::endl;
            iss >> tempColour;
            //set colour everytime encounter with 'usemtl'
            for (int i = 0; i < colourList.size(); ++i) {
                if (tempColour == colourList[i].name){
                    colour = colourList[i];
                }
            }
        }
    }

    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
//            std::cout << "position pixel: " << "(" << x << ", " << y << ")" << std::endl;
//            lookAt(glm::vec3(0, 0, 0));
            float focalLength = 2.0;
            uint32_t black = (0 << 24) + (0 << 16) + (0 << 8) + 0;
            uint32_t colour_t = black;
            glm::vec3 ray = computeRay(x, y);
            RayTriangleIntersection intersection;
            intersection = getClosestIntersection(cam, ray, coordinates);

            // get surface to light vector and its length
            glm::vec3 surfaceToLight = computeSurfaceToLight(intersection.intersectionPoint, lightSource);
            float length = surfaceToLightLength(surfaceToLight);

            // Calculate first reflection of light on surface for specular light
            glm::vec3 incidence = glm::normalize(-surfaceToLight);
            glm::vec3 normal = intersection.intersectedTriangle.normal;
            glm::vec3 reflect = incidence - (2 * normal) * (glm::dot(incidence, normal));
            // Dot product reflection and -camray
            float dp_vr = glm::dot(glm::normalize(reflect), -ray);
            dp_vr = pow(dp_vr, 64);

            // Calculate dp of normal and surface to light
            float dp_nl = glm::dot(intersection.intersectedTriangle.normal, glm::normalize(surfaceToLight));

//            colour_t = getColourShaded(length, dp_nl, colourList, intersection.intersectedTriangle);
//            window.setPixelColour(x, y, colour_t);

            // check for shadows
            bool shadow = checkShadow(intersection, intersection.intersectionPoint, glm::normalize(surfaceToLight), coordinates);
            if (shadow){
                // check if surface is facing camera
                float ang = glm::dot(intersection.intersectedTriangle.normal, ray);
                if (ang >= 0){
                    colour_t = getColourShaded(length, dp_nl, colourList, intersection.intersectedTriangle);
                    window.setPixelColour(x, y, colour_t);
                } else if (ang < 0) {
                    colour_t = getColourShadow(length, colourList, intersection.intersectedTriangle);
                    window.setPixelColour(x, y, colour_t);
                }
            } else if (!shadow) {
                float ang = glm::dot(intersection.intersectedTriangle.normal, glm::normalize(surfaceToLight));
                colour_t = getColourShaded(length, dp_nl, colourList, intersection.intersectedTriangle);
                window.setPixelColour(x, y, colour_t);
            }
//            }
        }
    }

}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_RIGHT) { //light right
            getLight(lightSource, glm::vec3(0.1, 0.0, 0.0));
//            depth.clear();
            window.clearPixels();
        } else if (event.key.keysym.sym == SDLK_LEFT) { //light left
            getLight(lightSource, glm::vec3(-0.1, 0.0, 0.0));
//            depth.clear();
            window.clearPixels();
        } else if (event.key.keysym.sym == SDLK_DOWN) { //light down
            getLight(lightSource, glm::vec3(0.0, -0.1, 0.0));
//            depth.clear();
            window.clearPixels();
        } else if (event.key.keysym.sym == SDLK_UP) { //light up
            getLight(lightSource, glm::vec3(0.0, 0.1, 0.0));
//            depth.clear();
            window.clearPixels();
        } else if (event.key.keysym.sym == 113) { // Q - cam forward
            getCamera(cam, glm::vec3(0.0, 0.0, -0.1));
            window.clearPixels();
        } else if (event.key.keysym.sym == 101) { // E - cam backward
            getCamera(cam, glm::vec3(0.0, 0.0, 0.1));
            window.clearPixels();
        } else if (event.key.keysym.sym == 114) { // R - Reset
            orbit = 0;
            if (mode == 0){
                cam = glm::vec3(0.0, 0.0, 4.0);
            } else if (mode == 1){
                cam = glm::vec3(0.0, 0.0, 8.0);
                lightSource = glm::vec3(0.0, 2.0, 0.0);
            } else if (mode == 2){
                cam = glm::vec3(0.0, 0.0, 4.0);
            }
            lookAt(glm::vec3(0, 0, 0));
            camOrientation = glm::mat3(-1, 0, 0,
                                       0, 1, 0,
                                       0, 0, 1);
            //            depth.clear();
            window.clearPixels();
        } else if (event.key.keysym.sym == 120) { // X - rotate about x-axis
            rotateCamera(cam, glm::mat3(1, 0, 0,
                                                    0, cos(0.017453), (-1 * sin(0.017453)),
                                                    0, sin(0.017453), cos(0.017453)));
            lookAt(glm::vec3(0, 0, 0));
            window.clearPixels();
        } else if (event.key.keysym.sym == 99) { // C - rotate about y-axis
            rotateCamera(cam, glm::mat3(cos(0.017453), 0, sin(0.017453),
                                                    0, 1, 0,
                                                    (-1 * sin(0.017453)), 0, cos(0.017453)));
            window.clearPixels();
        } else if (event.key.keysym.sym == 122) { // Z - rotate about z-axis
            rotateCamera(cam, glm::mat3(cos(0.017453), (-1 * sin(0.017453)), 0,
                                                    sin(0.017453),  cos(0.017453), 0,
                                                    0, 0, 1));
            window.clearPixels();
        } else if (event.key.keysym.sym == 116) { // T - rotate light about x-axis
            rotateLight(lightSource, glm::mat3(1, 0, 0,
                                        0, cos(0.017453), (-1 * sin(0.017453)),
                                        0, sin(0.017453), cos(0.017453)));
            lookAt(glm::vec3(0, 0, 0));
            window.clearPixels();
        } else if (event.key.keysym.sym == 121) { // Y - rotate light about y-axis
            rotateLight(lightSource, glm::mat3(cos(0.017453), 0, sin(0.017453),
                                        0, 1, 0,
                                        (-1 * sin(0.017453)), 0, cos(0.017453)));
            window.clearPixels();
        } else if (event.key.keysym.sym == 122) { // U - rotate light about z-axis
            rotateLight(lightSource, glm::mat3(cos(0.017453), (-1 * sin(0.017453)), 0,
                                        sin(0.017453),  cos(0.017453), 0,
                                        0, 0, 1));
            window.clearPixels();
        } else if (event.key.keysym.sym == 119) { // W camera up
            getCamera(cam, glm::vec3(0.0, 0.1, 0.0));
            window.clearPixels();
        } else if (event.key.keysym.sym == 97) { // A camera left
            getCamera(cam, glm::vec3(-0.1, 0.0, 0.0));
            window.clearPixels();
        } else if (event.key.keysym.sym == 115) { // S camera down
            getCamera(cam, glm::vec3(0.0, -0.1, 0.0));
            window.clearPixels();
        } else if (event.key.keysym.sym == 100) { // D camera right
            getCamera(cam, glm::vec3(0.1, 0.0, 0.0));
            window.clearPixels();
        } else if (event.key.keysym.sym == 49) { // 1 - camera zoom in
            getFocal(focalLength, 0.1);
            window.clearPixels();
        } else if (event.key.keysym.sym == 50) { // 2 - camera zoom out
            getFocal(focalLength, -0.1);
            window.clearPixels();
        } else if (event.key.keysym.sym == 108) { // L - lookAt()
            window.clearPixels();
            lookAt(glm::vec3(0, 0, 0));
        } else if (event.key.keysym.sym == 111) { // O - Orbit()
            window.clearPixels();
            if (orbit == 0){
                orbit = 1;
            } else if (orbit == 1){
                orbit = 0;
            }
        } else if (event.key.keysym.sym == 109) { // M - rendering mode
            window.clearPixels();
            if (mode == 0){
                getCamera(cam, glm::vec3(0.0, 0.0, 5.2));
                mode = 1;
                lookAt(glm::vec3(0, 0, 0));
            } else if (mode == 1){
                getCamera(cam, glm::vec3(0.0, 0.0, -5.2));
                mode = 2;
                lookAt(glm::vec3(0, 0, 0));
            } else if (mode == 2){
                mode = 0;
                lookAt(glm::vec3(0, 0, 0));
            }
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}


int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
    orbit = 0;
    mode = 1;
    int index = 0;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
        setDepthToZero();
//        window.clearPixels();
        if (orbit == 0){
            window.clearPixels();
            if (mode == 0){
                drawRaster(window,  "cornell-box.obj", "cornell-box.mtl", cam);
            } else if (mode == 1) {
                drawRayTrace(window,  "cornell-box.obj", "cornell-box.mtl", cam);
            } else if (mode == 2){
                drawWireframe(window,  "cornell-box.obj", "cornell-box.mtl", cam);
            }
        } else if (orbit == 1){
            window.clearPixels();
            rotateCamera(cam, glm::mat3(cos(0.017453), 0, sin(0.017453),
                                        0, 1, 0,
                                        (-1 * sin(0.017453)), 0, cos(0.017453)));
            lookAt(glm::vec3(0, 0, 0));
            if (mode == 0){
                drawRaster(window,  "cornell-box.obj", "cornell-box.mtl", cam);
            } else if (mode == 1) {
                drawRayTrace(window,  "cornell-box.obj", "cornell-box.mtl", cam);
            } else if (mode == 2){
                drawWireframe(window,  "cornell-box.obj", "cornell-box.mtl", cam);
            }
        }
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
//        window.savePPM("savedOutput/" + std::to_string(index) + ".ppm");
        index++;
	}

}
