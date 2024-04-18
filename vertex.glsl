#version 300 es
// Vertex coordinates in object space for the render quad
in vec3 vertexPosition;
uniform float canvasWidth;
uniform float canvasHeight;
uniform vec3 cameraPosition;
uniform mat3 cameraRotation;
uniform bool isOrthographicProjection;
uniform float orthographicFOV;
uniform float perspectiveFOV;

uniform mat4 mMatrix;
uniform mat4 vMatrix;
uniform mat4 pMatrix;


out vec3 origin;
out vec3 dir;

void main() {
    float aspectRatio = canvasWidth/canvasHeight;
    vec3 origin_camSpace, dir_camSpace;
    if (isOrthographicProjection) {
        origin_camSpace = vec3(vertexPosition.x*orthographicFOV*aspectRatio,
        vertexPosition.y*orthographicFOV,
        0);
        dir_camSpace = vec3(0, 0, -1);
    }
    else { // perspective projection
        //origin_camSpace = vec3(0);
        origin_camSpace =  mat3(vMatrix)*vec3(0.0, 0.0, 0.0);
        dir_camSpace = vec3(vertexPosition.x*aspectRatio,
        vertexPosition.y,
        -1.0/tan(radians(perspectiveFOV)));
    }
    origin = cameraPosition + cameraRotation*origin_camSpace;
    dir = normalize(cameraRotation*dir_camSpace);
    gl_Position = vec4(vertexPosition, 1.0);

}
