#version 300 es

precision highp float;

// A texture sampling unit, which is bound to the render quad texture buffer
// uniform sampler2D textureRendered;

// Texture coordinates coming from the vertex shader, interpolated through the rasterizer
// in vec2 fragmentTextureCoordinates;
out vec4 fragColor;

in vec3 origin;
in vec3 dir;

const int maxRayTraceDepth = 3;
uniform vec3 ambient;
uniform vec3 diffuse;
uniform vec3 specular;
uniform float shininess; // Specular shininess factor = q
uniform float lightIntensity;
uniform vec4 lightPosition;
uniform bool lightInCamspace;

const float k_fresnel = 0.4;

struct Sphere
{
    vec3 centre;
    float radius;
    vec3 colour;
};
struct Plane
{
    vec3 point;
    vec3 normal;
    vec3 colour;
};

struct Intersection
{
    bool flag; // false: not valid Intersection
    float u;
    vec3 position;
    vec3 normal;
    vec3 colour;
};

const int sphere_num = 6;
Sphere sphere[6];
Plane plane;

Intersection sphereIntersection(vec3 p0, vec3 d, Sphere thisSphere) {
    
    Intersection intersec;
    vec3 ps = thisSphere.centre;
    float r = thisSphere.radius;
    vec3 delta_p = p0 - ps;
    float temp = pow(dot(d, delta_p), 2.0) - pow(length(delta_p), 2.0) + r*r;
    float u1 = -dot(d, delta_p) - sqrt(temp);
    float u2 = -dot(d, delta_p) + sqrt(temp);


    if((temp <= 0.001 || u1 <= 0.001)) {
        intersec.flag = false;
        return intersec;
    }

    // calculate the Intersection
    intersec.flag = true;
    intersec.u = u1;
    intersec.position = p0 + u1 * d;
    intersec.normal = normalize(intersec.position - ps);
    intersec.colour = thisSphere.colour;
    return intersec;
}

Intersection planeIntersection(vec3 p0, vec3 d, Plane thisPlane) {
    Intersection intersec;
    vec3 p1 = thisPlane.point;
    vec3 n = thisPlane.normal;
    float u = -dot(p0-p1, n)/dot(d, n);
    if(u <= 0.001) { // add offset to avoid self-shadowing
        intersec.flag = false;
        return intersec;
    }

    // calculate the Intersection
    intersec.flag = true;
    intersec.u = u;
    intersec.position = p0 + u * d;
    intersec.normal = n;
    float modbase = 1.0;
    float modx = mod(intersec.position.x, modbase);
    float modz = mod(intersec.position.z, modbase);

    // checkerboard pattern results even-odd cross grid
    if((modx<modbase/2.0 && modz<modbase/2.0) || (modx>=modbase/2.0 && modz>=modbase/2.0))
        intersec.colour = vec3(0.7, 0.7, 0.7);
    else
        intersec.colour = vec3(0.3, 0.3, 0.3);

    return intersec;
}

// Shadow
bool shadow_check(Intersection intersection) {
    // slightly move the ray origin outwards of the object along the surface normal
    vec3 p0 = vec3(lightPosition);
    vec3 d = normalize(intersection.position - vec3(lightPosition));

    float u_this = length(intersection.position - vec3(lightPosition));

    float min_u = 100000.0;
    bool flag_sphere = false;
    for(int i=0; i<sphere_num; i++) {
        Intersection intersectionSphere = sphereIntersection(p0, d, sphere[i]);
        if(intersectionSphere.flag == true) {
            flag_sphere = true;
            min_u = min(min_u, intersectionSphere.u);
            //break;
        }
    }
    Intersection intersectionPlane = planeIntersection(p0, d, plane);
    if(intersectionPlane.flag == true) {
        min_u = min(min_u, intersectionPlane.u);
    }

    if(intersectionPlane.flag == false && flag_sphere == false) {
        return false;
    }

    if(min_u < u_this-0.001) {
        return true;
    }
    return false;
}

// light ads
vec3 ads(Intersection intersection)
{
    vec3 outPosition = intersection.position;
	vec3 n = intersection.normal;
	vec3 l;
    float d;
	if(lightInCamspace == true) {
        l = normalize(-outPosition);
        d = length(outPosition);
    }
	else {
        l = normalize(vec3(lightPosition) - outPosition);
        d = length(vec3(lightPosition) - outPosition);
    }
		
	vec3 v = normalize(-outPosition);
	vec3 r = reflect(-l, n);

    float new_i = lightIntensity/(4.0*3.1415926*d*d);

	vec3 res =  ambient + 
                new_i * (diffuse * max(dot(l, n), 0.0) + 
                specular * pow(max(dot(r,v), 0.0), shininess));
    return vec3(res);
}


struct Ray
{
    vec3 p0;
    vec3 d;
    float intensity_k;
};

// Refraction - Reflection
// sphere: 0.4-0.6, plane: 0.0-1.0

Ray Refract(Ray ray, Intersection intersection) {
    Ray new_ray;
    new_ray.p0 = intersection.position - intersection.normal*0.002; // another side
    new_ray.d = normalize(refract(ray.d, intersection.normal, 0.6));
    new_ray.intensity_k = ray.intensity_k*(1.0-k_fresnel);
    return new_ray;
}

Ray Reflect(Ray ray, Intersection intersection) {
    Ray new_ray;
    new_ray.p0 = intersection.position + intersection.normal*0.001;
    new_ray.d = normalize(reflect(ray.d, intersection.normal));
    new_ray.intensity_k = ray.intensity_k * k_fresnel;
    return new_ray;
}

Ray defalult_ray() {
    Ray rayx;
    rayx.p0 = vec3(0.0);
    rayx.d = vec3(0.0);
    rayx.intensity_k = 0.0;
    return rayx;
}

//const int MAX_LIST_SIZE = 16;
Ray ray_list[4];
vec3 Result[7];
int res_pos = 0;


// Main program for each fragment of the render quad
void main() {
    sphere[0].centre = vec3(-2.0, 1.5, -3.5);
    sphere[0].radius = 1.5;
    sphere[0].colour = vec3(0.8,0.8,0.8);
    sphere[1].centre = vec3(-0.5, 0.0, -2.0);
    sphere[1].radius = 0.6;
    sphere[1].colour = vec3(0.3,0.8,0.3);
    sphere[2].centre = vec3(1.0, 0.7, -2.2);
    sphere[2].radius = 0.8;
    sphere[2].colour = vec3(0.3,0.8,0.8);
    sphere[3].centre = vec3(0.7, -0.3, -1.2);
    sphere[3].radius = 0.2;
    sphere[3].colour = vec3(0.8,0.8,0.3);
    sphere[4].centre = vec3(-0.7, -0.3, -1.2);
    sphere[4].radius = 0.2;
    sphere[4].colour = vec3(0.8,0.3,0.3);
    sphere[5].centre = vec3(0.2, -0.2, -1.2);
    sphere[5].radius = 0.3;
    sphere[5].colour = vec3(0.8,0.3,0.8);
    plane.point = vec3(0, -0.5, 0);
    plane.normal = vec3(0, 1.0, 0);
    plane.colour = vec3(1, 1, 1);
    //scene definition end

    vec3 RayTracecolour = vec3(0.0, 0.0, 0.0);
    int now_list_size = 1;
    ray_list[0] = Ray(origin, normalize(dir), 1.0);

    for(int depth = 0; depth < maxRayTraceDepth; depth++) {
        int last_list_size = now_list_size;
        now_list_size = 0;
        Ray new_ray_list[4];

        for(int it_pos=0;it_pos<last_list_size;it_pos++) {

            Ray ray = ray_list[it_pos];
            vec3 p0 = ray_list[it_pos].p0;
            vec3 d = ray_list[it_pos].d;
            if(ray.intensity_k < 0.001) {
                Result[res_pos] = vec3(0.0, 0.0, 0.0);
                res_pos++;
                if(depth != maxRayTraceDepth-1) {
                    new_ray_list[now_list_size] = defalult_ray();
                    new_ray_list[now_list_size+1] = defalult_ray();
                    now_list_size += 2;
                }

                continue;
            }

            /*----------------------------------*/

            Intersection intersectionList[7];

            for(int i=0; i<sphere_num; i++) {
                intersectionList[i] = sphereIntersection(p0, d, sphere[i]);
            }
            intersectionList[sphere_num] = planeIntersection(p0, d, plane);
            Intersection closeIntersection;
            float min_u;
            int valid_flag = -1;
            for(int i=0; i<sphere_num+1; i++) {
                if(intersectionList[i].flag == true) {
                    if(valid_flag == -1) {
                        closeIntersection = intersectionList[i];
                        min_u = intersectionList[i].u;
                        valid_flag = i;
                    }
                    else if(intersectionList[i].u < min_u) {
                        closeIntersection = intersectionList[i];
                        min_u = intersectionList[i].u;
                        valid_flag = i;
                    }
                    else {
                        ;
                    }
                }
            }

            if(valid_flag == -1) { // no valid intersection
                    Result[res_pos] = vec3(0.0);
                    res_pos++;
                    if(depth != maxRayTraceDepth-1) {
                        new_ray_list[now_list_size] = defalult_ray();
                        new_ray_list[now_list_size+1] = defalult_ray();
                        now_list_size += 2;
                    }
                    continue;
            }

            vec3 new_add_color = vec3(0.0, 0.0, 0.0);
            if(valid_flag == sphere_num && shadow_check(closeIntersection) == true) {
                Result[res_pos] = vec3(ambient);
                res_pos++;
                if(depth != maxRayTraceDepth-1) {
                    new_ray_list[now_list_size] = defalult_ray();
                    new_ray_list[now_list_size+1] = defalult_ray();
                    now_list_size += 2;
                }
            }
            else {
                if(valid_flag == sphere_num) { // on plane, only reflect
                    new_add_color += ads(closeIntersection) * closeIntersection.colour;
                    Ray new_ray_reflect = Reflect(ray, closeIntersection);
                    Result[res_pos] = new_add_color;
                    res_pos++;
                    if(depth != maxRayTraceDepth-1) {
                        new_ray_list[now_list_size] = new_ray_reflect;
                        new_ray_list[now_list_size+1] = defalult_ray();
                        now_list_size += 2;
                    }
                }
                else { // on sphere
                    Ray new_ray_reflect = Reflect(ray, closeIntersection);
                    Ray new_ray_refract = Refract(ray, closeIntersection);
                    new_add_color = ads(closeIntersection) * closeIntersection.colour;
                    Result[res_pos] = new_add_color;
                    res_pos++;
                    if(depth != maxRayTraceDepth-1) {
                        new_ray_list[now_list_size] = new_ray_reflect;
                        new_ray_list[now_list_size+1] = new_ray_refract;
                        now_list_size += 2;
                    }
                }
                
            }
        }
        ray_list = new_ray_list;
        //now_list_size = min(MAX_LIST_SIZE, now_list_size);
    }

    for(int i=maxRayTraceDepth-1;i>0;i--) {
        int start = (1<<i)-1;
        for(int j=start;j<start*2;j+=2) {
            Result[(j-1)/2] += (Result[j] + Result[j+1]); 
        }
    }
    RayTracecolour = Result[0];
    fragColor = vec4(RayTracecolour, 1.0);
}

