

// ---------------------------------------------------------------------- //
// External code
// ---------------------------------------------------------------------- //
// http://iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdCylinder(vec3 p, vec3 c) {return length(p.xz-c.xy)-c.z;}
float sdSphere(vec3 p, float s) {return length(p)-s;}
// http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/
mat4 rotationMatrix(vec3 axis, float angle) {
    axis = normalize(axis);float s = sin(angle);float c = cos(angle);float oc = 1.0 - c;
    return mat4(oc*axis.x*axis.x+c,       oc*axis.x*axis.y-axis.z*s,oc*axis.z*axis.x+axis.y*s,0.0,
                oc*axis.x*axis.y+axis.z*s,oc*axis.y*axis.y+c,       oc*axis.y*axis.z-axis.x*s,0.0,
                oc*axis.z*axis.x-axis.y*s,oc*axis.y*axis.z+axis.x*s,oc*axis.z*axis.z+c,       0.0,
                0.0,0.0,0.0,1.0);
}

////////////////////////////////////////////////////////////////
//
//                           HG_SDF
//
//     GLSL LIBRARY FOR BUILDING SIGNED DISTANCE BOUNDS
//
//     version 2016-01-04
//
//     Check http://mercury.sexy/hg_sdf for updates
//     and usage examples. Send feedback to spheretracing@mercury.sexy.
//
//     Brought to you by MERCURY http://mercury.sexy
//
//
//
// Released as Creative Commons Attribution-NonCommercial (CC BY-NC)
//
////////////////////////////////////////////////////////////////
//
// How to use this:
//
// 1. Build some system to #include glsl files in each other.
//   Include this one at the very start. Or just paste everywhere.
// 2. Build a sphere tracer. See those papers:
//   * "Sphere Tracing" http://graphics.cs.illinois.edu/sites/default/files/zeno.pdf
//   * "Enhanced Sphere Tracing" http://lgdv.cs.fau.de/get/2234
//   The Raymnarching Toolbox Thread on pouet can be helpful as well
//   http://www.pouet.net/topic.php?which=7931&page=1
//   and contains links to many more resources.
// 3. Use the tools in this library to build your distance bound f().
// 4. ???
// 5. Win a compo.
// 
// (6. Buy us a beer or a good vodka or something, if you like.)
//
////////////////////////////////////////////////////////////////
//
// Table of Contents:
//
// * Helper functions and macros
// * Collection of some primitive objects
// * Domain Manipulation operators
// * Object combination operators
//
////////////////////////////////////////////////////////////////
//
// Why use this?
//
// The point of this lib is that everything is structured according
// to patterns that we ended up using when building geometry.
// It makes it more easy to write code that is reusable and that somebody
// else can actually understand. Especially code on Shadertoy (which seems
// to be what everybody else is looking at for "inspiration") tends to be
// really ugly. So we were forced to do something about the situation and
// release this lib ;)
//
// Everything in here can probably be done in some better way.
// Please experiment. We'd love some feedback, especially if you
// use it in a scene production.
//
// The main patterns for building geometry this way are:
// * Stay Lipschitz continuous. That means: don't have any distance
//   gradient larger than 1. Try to be as close to 1 as possible -
//   Distances are euclidean distances, don't fudge around.
//   Underestimating distances will happen. That's why calling
//   it a "distance bound" is more correct. Don't ever multiply
//   distances by some value to "fix" a Lipschitz continuity
//   violation. The invariant is: each fSomething() function returns
//   a correct distance bound.
// * Use very few primitives and combine them as building blocks
//   using combine opertors that preserve the invariant.
// * Multiply objects by repeating the domain (space).
//   If you are using a loop inside your distance function, you are
//   probably doing it wrong (or you are building boring fractals).
// * At right-angle intersections between objects, build a new local
//   coordinate system from the two distances to combine them in
//   interesting ways.
// * As usual, there are always times when it is best to not follow
//   specific patterns.
//
////////////////////////////////////////////////////////////////
//
// FAQ
//
// Q: Why is there no sphere tracing code in this lib?
// A: Because our system is way too complex and always changing.
//    This is the constant part. Also we'd like everyone to
//    explore for themselves.
//
// Q: This does not work when I paste it into Shadertoy!!!!
// A: Yes. It is GLSL, not GLSL ES. We like real OpenGL
//    because it has way more features and is more likely
//    to work compared to browser-based WebGL. We recommend
//    you consider using OpenGL for your productions. Most
//    of this can be ported easily though.
//
// Q: How do I material?
// A: We recommend something like this:
//    Write a material ID, the distance and the local coordinate
//    p into some global variables whenever an object's distance is
//    smaller than the stored distance. Then, at the end, evaluate
//    the material to get color, roughness, etc., and do the shading.
//
// Q: I found an error. Or I made some function that would fit in
//    in this lib. Or I have some suggestion.
// A: Awesome! Drop us a mail at spheretracing@mercury.sexy.
//
// Q: Why is this not on github?
// A: Because we were too lazy. If we get bugged about it enough,
//    we'll do it.
//
// Q: Your license sucks for me.
// A: Oh. What should we change it to?
//
// Q: I have trouble understanding what is going on with my distances.
// A: Some visualization of the distance field helps. Try drawing a
//    plane that you can sweep through your scene with some color
//    representation of the distance field at each point and/or iso
//    lines at regular intervals. Visualizing the length of the
//    gradient (or better: how much it deviates from being equal to 1)
//    is immensely helpful for understanding which parts of the
//    distance field are broken.
//
////////////////////////////////////////////////////////////////






////////////////////////////////////////////////////////////////
//
//             HELPER FUNCTIONS/MACROS
//
////////////////////////////////////////////////////////////////

#define PI 3.14159265
#define TAU (2*PI)
#define PHI (sqrt(5.)*0.5 + 0.5)

// Clamp to [0,1] - this operation is free under certain circumstances.
// For further information see
// http://www.humus.name/Articles/Persson_LowLevelThinking.pdf and
// http://www.humus.name/Articles/Persson_LowlevelShaderOptimization.pdf
#define saturate(x) clamp(x, 0., 1.)

// Sign function that doesn't return 0
float sgn(float x) {
    return (x<0.)?-1.:1.;
}

vec2 sgn(vec2 x) {
    return vec2(sgn(x.x), sgn(x.y));
}

float square (float x) {
    return x*x;
}

vec2 square (vec2 x) {
    return x*x;
}

vec3 square (vec3 x) {
    return x*x;
}

float lengthSqr(vec3 x) {
    return dot(x, x);
}


// Maximum/minumum elements of a vector
float vmax(vec2 v) {
    return max(v.x, v.y);
}

float vmax(vec3 v) {
    return max(max(v.x, v.y), v.z);
}

float vmax(vec4 v) {
    return max(max(v.x, v.y), max(v.z, v.w));
}

float vmin(vec2 v) {
    return min(v.x, v.y);
}

float vmin(vec3 v) {
    return min(min(v.x, v.y), v.z);
}

float vmin(vec4 v) {
    return min(min(v.x, v.y), min(v.z, v.w));
}




////////////////////////////////////////////////////////////////
//
//             PRIMITIVE DISTANCE FUNCTIONS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that is a distance function is called fSomething.
// The first argument is always a point in 2 or 3-space called <p>.
// Unless otherwise noted, (if the object has an intrinsic "up"
// side or direction) the y axis is "up" and the object is
// centered at the origin.
//
////////////////////////////////////////////////////////////////

float fSphere(vec3 p, float r) {
    return length(p) - r;
}

// Plane with normal n (n is normalized) at some distance from the origin
float fPlane(vec3 p, vec3 n, float distanceFromOrigin) {
    return dot(p, n) + distanceFromOrigin;
}

// Cheap Box: distance to corners is overestimated
float fBoxCheap(vec3 p, vec3 b) { //cheap box
    return vmax(abs(p) - b);
}

// Box: correct distance to corners
float fBox(vec3 p, vec3 b) {
    vec3 d = abs(p) - b;
    return length(max(d, vec3(0))) + vmax(min(d, vec3(0)));
}

// Same as above, but in two dimensions (an endless box)
float fBox2Cheap(vec2 p, vec2 b) {
    return vmax(abs(p)-b);
}

float fBox2(vec2 p, vec2 b) {
    vec2 d = abs(p) - b;
    return length(max(d, vec2(0))) + vmax(min(d, vec2(0)));
}


// Endless "corner"
float fCorner (vec2 p) {
    return length(max(p, vec2(0))) + vmax(min(p, vec2(0)));
}

// Blobby ball object. You've probably seen it somewhere. This is not a correct distance bound, beware.
float fBlob(vec3 p) {
    p = abs(p);
    if (p.x < max(p.y, p.z)) p = p.yzx;
    if (p.x < max(p.y, p.z)) p = p.yzx;
    float b = max(max(max(
        dot(p, normalize(vec3(1, 1, 1))),
        dot(p.xz, normalize(vec2(PHI+1., 1.)))),
        dot(p.yx, normalize(vec2(1., PHI)))),
        dot(p.xz, normalize(vec2(1., PHI))));
    float l = length(p);
    return l - 1.5 - 0.2 * (1.5 / 2.)* cos(min(sqrt(1.01 - b / l)*(PI / 0.25), PI));
}

// Cylinder standing upright on the xz plane
float fCylinder(vec3 p, float r, float height) {
    float d = length(p.xz) - r;
    d = max(d, abs(p.y) - height);
    return d;
}

// Capsule: A Cylinder with round caps on both sides
float fCapsule(vec3 p, float r, float c) {
    return mix(length(p.xz) - r, length(vec3(p.x, abs(p.y) - c, p.z)) - r, step(c, abs(p.y)));
}

// Distance to line segment between <a> and <b>, used for fCapsule() version 2below
float fLineSegment(vec3 p, vec3 a, vec3 b) {
    vec3 ab = b - a;
    float t = saturate(dot(p - a, ab) / dot(ab, ab));
    return length((ab*t + a) - p);
}

// Capsule version 2: between two end points <a> and <b> with radius r 
float fCapsule(vec3 p, vec3 a, vec3 b, float r) {
    return fLineSegment(p, a, b) - r;
}

// Torus in the XZ-plane
float fTorus(vec3 p, float smallRadius, float largeRadius) {
    return length(vec2(length(p.xz) - largeRadius, p.y)) - smallRadius;
}

// A circle line. Can also be used to make a torus by subtracting the smaller radius of the torus.
float fCircle(vec3 p, float r) {
    float l = length(p.xz) - r;
    return length(vec2(p.y, l));
}

// A circular disc with no thickness (i.e. a cylinder with no height).
// Subtract some value to make a flat disc with rounded edge.
float fDisc(vec3 p, float r) {
    float l = length(p.xz) - r;
    return l < 0. ? abs(p.y) : length(vec2(p.y, l));
}

// Hexagonal prism, circumcircle variant
float fHexagonCircumcircle(vec3 p, vec2 h) {
    vec3 q = abs(p);
    return max(q.y - h.y, max(q.x*sqrt(3.)*0.5 + q.z*0.5, q.z) - h.x);
    //this is mathematically equivalent to this line, but less efficient:
    //return max(q.y - h.y, max(dot(vec2(cos(PI/3), sin(PI/3)), q.zx), q.z) - h.x);
}

// Hexagonal prism, incircle variant
float fHexagonIncircle(vec3 p, vec2 h) {
    return fHexagonCircumcircle(p, vec2(h.x*sqrt(3.)*0.5, h.y));
}

// Cone with correct distances to tip and base circle. Y is up, 0 is in the middle of the base.
float fCone(vec3 p, float radius, float height) {
    vec2 q = vec2(length(p.xz), p.y);
    vec2 tip = q - vec2(0, height);
    vec2 mantleDir = normalize(vec2(height, radius));
    float mantle = dot(tip, mantleDir);
    float d = max(mantle, -q.y);
    float projected = dot(tip, vec2(mantleDir.y, -mantleDir.x));
    
    // distance to tip
    if ((q.y > height) && (projected < 0.)) {
        d = max(d, length(tip));
    }
    
    // distance to base ring
    if ((q.x > radius) && (projected > length(vec2(height, radius)))) {
        d = max(d, length(q - vec2(radius, 0)));
    }
    return d;
}

//
// "Generalized Distance Functions" by Akleman and Chen.
// see the Paper at https://www.viz.tamu.edu/faculty/ergun/research/implicitmodeling/papers/sm99.pdf
//
// This set of constants is used to construct a large variety of geometric primitives.
// Indices are shifted by 1 compared to the paper because we start counting at Zero.
// Some of those are slow whenever a driver decides to not unroll the loop,
// which seems to happen for fIcosahedron und fTruncatedIcosahedron on nvidia 350.12 at least.
// Specialized implementations can well be faster in all cases.
//

vec3 GDFVectors[19];

void hg_sdf_init() {
    GDFVectors[0] = normalize(vec3(1, 0, 0));
    GDFVectors[1] = normalize(vec3(0, 1, 0));
    GDFVectors[2] = normalize(vec3(0, 0, 1));

    GDFVectors[3] = normalize(vec3(1, 1, 1 ));
    GDFVectors[4] = normalize(vec3(-1, 1, 1));
    GDFVectors[5] = normalize(vec3(1, -1, 1));
    GDFVectors[6] = normalize(vec3(1, 1, -1));

    GDFVectors[7] = normalize(vec3(0., 1., PHI+1.));
    GDFVectors[8] = normalize(vec3(0., -1., PHI+1.));
    GDFVectors[9] = normalize(vec3(PHI+1., 0., 1.));
    GDFVectors[10] = normalize(vec3(-PHI-1., 0., 1.));
    GDFVectors[11] = normalize(vec3(1., PHI+1., 0.));
    GDFVectors[12] = normalize(vec3(-1., PHI+1., 0.));

    GDFVectors[13] = normalize(vec3(0, PHI, 1));
    GDFVectors[14] = normalize(vec3(0, -PHI, 1));
    GDFVectors[15] = normalize(vec3(1, 0, PHI));
    GDFVectors[16] = normalize(vec3(-1, 0, PHI));
    GDFVectors[17] = normalize(vec3(PHI, 1, 0));
    GDFVectors[18] = normalize(vec3(-PHI, 1, 0));
}
//);

// Version with variable exponent.
// This is slow and does not produce correct distances, but allows for bulging of objects.
#define fGDF_t(begin,end) float d = 0.; for (int i = begin; i <= end; ++i) {d += pow(abs(dot(p, GDFVectors[i])), e);} return pow(d, 1./e) - r;
// Original method, not suitable for WebGL:
float fGDF(vec3 p, float r, float e, int begin, int end) {
    float d = 0.;
    for (int i = begin; i <= end; ++i)
        d += pow(abs(dot(p, GDFVectors[i])), e);
    return pow(d, 1./e) - r;
}
float fGDF36(vec3 p, float r, float e) { fGDF_t(3,6) }
float fGDF1318(vec3 p, float r, float e) { fGDF_t(13,18) }
float fGDF318(vec3 p, float r, float e) { fGDF_t(3,18) }
float fGDF312(vec3 p, float r, float e) { fGDF_t(3,12) }
float fGDF06(vec3 p, float r, float e) { fGDF_t(0,6) }
#undef fGDF_t

// Version with without exponent, creates objects with sharp edges and flat faces
#define fGDF_t(begin,end) float d = 0.; for (int i = begin; i <= end; ++i){ d = max(d, abs(dot(p, GDFVectors[i])));} return d - r;
// Original method, not suitable for WebGL:
float fGDF(vec3 p, float r, int begin, int end) {
    float d = 0.;
    for (int i = begin; i <= end; ++i)
        d = max(d, abs(dot(p, GDFVectors[i])));
    return d - r;
}
float fGDF36(vec3 p, float r) { fGDF_t(3,6) }
float fGDF1318(vec3 p, float r) { fGDF_t(13,18) }
float fGDF318(vec3 p, float r) { fGDF_t(3,18) }
float fGDF312(vec3 p, float r) { fGDF_t(3,12) }
float fGDF06(vec3 p, float r) { fGDF_t(0,6) }
#undef fGDF_t

// Primitives follow:

float fOctahedron(vec3 p, float r, float e) {
    return fGDF36(p, r, e);
}

float fDodecahedron(vec3 p, float r, float e) {
    return fGDF1318(p, r, e);
}

float fIcosahedron(vec3 p, float r, float e) {
    return fGDF312(p, r, e);
}

float fTruncatedOctahedron(vec3 p, float r, float e) {
    return fGDF06(p, r, e);
}

float fTruncatedIcosahedron(vec3 p, float r, float e) {
    return fGDF318(p, r, e);
}

float fOctahedron(vec3 p, float r) {
    return fGDF36(p, r);
}

float fDodecahedron(vec3 p, float r) {
    return fGDF1318(p, r);
}

float fIcosahedron(vec3 p, float r) {
    return fGDF312(p, r);
}

float fTruncatedOctahedron(vec3 p, float r) {
    return fGDF06(p, r);
}

float fTruncatedIcosahedron(vec3 p, float r) {
    return fGDF318(p, r);
}


////////////////////////////////////////////////////////////////
//
//                DOMAIN MANIPULATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that modifies the domain is named pSomething.
//
// Many operate only on a subset of the three dimensions. For those,
// you must choose the dimensions that you want manipulated
// by supplying e.g. <p.x> or <p.zx>
//
// <inout p> is always the first argument and modified in place.
//
// Many of the operators partition space into cells. An identifier
// or cell index is returned, if possible. This return value is
// intended to be optionally used e.g. as a random seed to change
// parameters of the distance functions inside the cells.
//
// Unless stated otherwise, for cell index 0, <p> is unchanged and cells
// are centered on the origin so objects don't have to be moved to fit.
//
//
////////////////////////////////////////////////////////////////



// Rotate around a coordinate axis (i.e. in a plane perpendicular to that axis) by angle <a>.
// Read like this: R(p.xz, a) rotates "x towards z".
// This is fast if <a> is a compile-time constant and slower (but still practical) if not.
void pR(inout vec2 p, float a) {
    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// Shortcut for 45-degrees rotation
void pR45(inout vec2 p) {
    p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

// Repeat space along one axis. Use like this to repeat along the x axis:
// <float cell = pMod1(p.x,5);> - using the return value is optional.
float pMod1(inout float p, float size) {
    float halfsize = size*0.5;
    float c = floor((p + halfsize)/size);
    p = mod(p + halfsize, size) - halfsize;
    return c;
}

// Same, but mirror every second cell so they match at the boundaries
float pModMirror1(inout float p, float size) {
    float halfsize = size*0.5;
    float c = floor((p + halfsize)/size);
    p = mod(p + halfsize,size) - halfsize;
    p *= mod(c, 2.0)*2. - 1.;
    return c;
}

// Repeat the domain only in positive direction. Everything in the negative half-space is unchanged.
float pModSingle1(inout float p, float size) {
    float halfsize = size*0.5;
    float c = floor((p + halfsize)/size);
    if (p >= 0.)
        p = mod(p + halfsize, size) - halfsize;
    return c;
}

// Repeat only a few times: from indices <start> to <stop> (similar to above, but more flexible)
float pModInterval1(inout float p, float size, float start, float stop) {
    float halfsize = size*0.5;
    float c = floor((p + halfsize)/size);
    p = mod(p+halfsize, size) - halfsize;
    if (c > stop) { //yes, this might not be the best thing numerically.
        p += size*(c - stop);
        c = stop;
    }
    if (c <start) {
        p += size*(c - start);
        c = start;
    }
    return c;
}


// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolar(inout vec2 p, float repetitions) {
    float angle = 2.*PI/repetitions;
    float a = atan(p.y, p.x) + angle/2.;
    float r = length(p);
    float c = floor(a/angle);
    a = mod(a,angle) - angle/2.;
    p = vec2(cos(a), sin(a))*r;
    // For an odd number of repetitions, fix cell index of the cell in -x direction
    // (cell index would be e.g. -5 and 5 in the two halves of the cell):
    if (abs(c) >= (repetitions/2.)) c = abs(c);
    return c;
}

// Repeat in two dimensions
vec2 pMod2(inout vec2 p, vec2 size) {
    vec2 c = floor((p + size*0.5)/size);
    p = mod(p + size*0.5,size) - size*0.5;
    return c;
}

// Same, but mirror every second cell so all boundaries match
vec2 pModMirror2(inout vec2 p, vec2 size) {
    vec2 halfsize = size*0.5;
    vec2 c = floor((p + halfsize)/size);
    p = mod(p + halfsize, size) - halfsize;
    p *= mod(c,vec2(2))*2. - vec2(1);
    return c;
}

// Same, but mirror every second cell at the diagonal as well
vec2 pModGrid2(inout vec2 p, vec2 size) {
    vec2 c = floor((p + size*0.5)/size);
    p = mod(p + size*0.5, size) - size*0.5;
    p *= mod(c,vec2(2))*2. - vec2(1);
    p -= size/2.;
    if (p.x > p.y) p.xy = p.yx;
    return floor(c/2.);
}

// Repeat in three dimensions
vec3 pMod3(inout vec3 p, vec3 size) {
    vec3 c = floor((p + size*0.5)/size);
    p = mod(p + size*0.5, size) - size*0.5;
    return c;
}

// Mirror at an axis-aligned plane which is at a specified distance <dist> from the origin.
float pMirror (inout float p, float dist) {
    float s = sgn(p);
    p = abs(p)-dist;
    return s;
}

// Mirror in both dimensions and at the diagonal, yielding one eighth of the space.
// translate by dist before mirroring.
vec2 pMirrorOctant (inout vec2 p, vec2 dist) {
    vec2 s = sgn(p);
    pMirror(p.x, dist.x);
    pMirror(p.y, dist.y);
    if (p.y > p.x)
        p.xy = p.yx;
    return s;
}

// Reflect space at a plane
float pReflect(inout vec3 p, vec3 planeNormal, float offset) {
    float t = dot(p, planeNormal)+offset;
    if (t < 0.) {
        p = p - (2.*t)*planeNormal;
    }
    return sgn(t);
}


////////////////////////////////////////////////////////////////
//
//             OBJECT COMBINATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// We usually need the following boolean operators to combine two objects:
// Union: OR(a,b)
// Intersection: AND(a,b)
// Difference: AND(a,!b)
// (a and b being the distances to the objects).
//
// The trivial implementations are min(a,b) for union, max(a,b) for intersection
// and max(a,-b) for difference. To combine objects in more interesting ways to
// produce rounded edges, chamfers, stairs, etc. instead of plain sharp edges we
// can use combination operators. It is common to use some kind of "smooth minimum"
// instead of min(), but we don't like that because it does not preserve Lipschitz
// continuity in many cases.
//
// Naming convention: since they return a distance, they are called fOpSomething.
// The different flavours usually implement all the boolean operators above
// and are called fOpUnionRound, fOpIntersectionRound, etc.
//
// The basic idea: Assume the object surfaces intersect at a right angle. The two
// distances <a> and <b> constitute a new local two-dimensional coordinate system
// with the actual intersection as the origin. In this coordinate system, we can
// evaluate any 2D distance function we want in order to shape the edge.
//
// The operators below are just those that we found useful or interesting and should
// be seen as examples. There are infinitely more possible operators.
//
// They are designed to actually produce correct distances or distance bounds, unlike
// popular "smooth minimum" operators, on the condition that the gradients of the two
// SDFs are at right angles. When they are off by more than 30 degrees or so, the
// Lipschitz condition will no longer hold (i.e. you might get artifacts). The worst
// case is parallel surfaces that are close to each other.
//
// Most have a float argument <r> to specify the radius of the feature they represent.
// This should be much smaller than the object size.
//
// Some of them have checks like "if ((-a < r) && (-b < r))" that restrict
// their influence (and computation cost) to a certain area. You might
// want to lift that restriction or enforce it. We have left it as comments
// in some cases.
//
// usage example:
//
// float fTwoBoxes(vec3 p) {
//   float box0 = fBox(p, vec3(1));
//   float box1 = fBox(p-vec3(1), vec3(1));
//   return fOpUnionChamfer(box0, box1, 0.2);
// }
//
////////////////////////////////////////////////////////////////


// The "Chamfer" flavour makes a 45-degree chamfered edge (the diagonal of a square of size <r>):
float fOpUnionChamfer(float a, float b, float r) {
    return min(min(a, b), (a - r + b)*sqrt(0.5));
}

// Intersection has to deal with what is normally the inside of the resulting object
// when using union, which we normally don't care about too much. Thus, intersection
// implementations sometimes differ from union implementations.
float fOpIntersectionChamfer(float a, float b, float r) {
    return max(max(a, b), (a + r + b)*sqrt(0.5));
}

// Difference can be built from Intersection or Union:
float fOpDifferenceChamfer (float a, float b, float r) {
    return fOpIntersectionChamfer(a, -b, r);
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float fOpUnionRound(float a, float b, float r) {
    vec2 u = max(vec2(r - a,r - b), vec2(0));
    return max(r, min (a, b)) - length(u);
}

float fOpIntersectionRound(float a, float b, float r) {
    vec2 u = max(vec2(r + a,r + b), vec2(0));
    return min(-r, max (a, b)) + length(u);
}

float fOpDifferenceRound (float a, float b, float r) {
    return fOpIntersectionRound(a, -b, r);
}


// The "Columns" flavour makes n-1 circular columns at a 45 degree angle:
float fOpUnionColumns(float a, float b, float r, float n) {
    if ((a < r) && (b < r)) {
        vec2 p = vec2(a, b);
        float columnradius = r*sqrt(2.)/((n-1.)*2.+sqrt(2.));
        pR45(p);
        p.x -= sqrt(2.)/2.*r;
        p.x += columnradius*sqrt(2.);
        if (mod(n,2.) == 1.) {
            p.y += columnradius;
        }
        // At this point, we have turned 45 degrees and moved at a point on the
        // diagonal that we want to place the columns on.
        // Now, repeat the domain along this direction and place a circle.
        pMod1(p.y, columnradius*2.);
        float result = length(p) - columnradius;
        result = min(result, p.x);
        result = min(result, a);
        return min(result, b);
    } else {
        return min(a, b);
    }
}

float fOpDifferenceColumns(float a, float b, float r, float n) {
    a = -a;
    float m = min(a, b);
    //avoid the expensive computation where not needed (produces discontinuity though)
    if ((a < r) && (b < r)) {
        vec2 p = vec2(a, b);
        float columnradius = r*sqrt(2.)/n/2.0;
        columnradius = r*sqrt(2.)/((n-1.)*2.+sqrt(2.));

        pR45(p);
        p.y += columnradius;
        p.x -= sqrt(2.)/2.*r;
        p.x += -columnradius*sqrt(2.)/2.;

        if (mod(n,2.) == 1.) {
            p.y += columnradius;
        }
        pMod1(p.y,columnradius*2.);

        float result = -length(p) + columnradius;
        result = max(result, p.x);
        result = min(result, a);
        return -min(result, b);
    } else {
        return -m;
    }
}

float fOpIntersectionColumns(float a, float b, float r, float n) {
    return fOpDifferenceColumns(a,-b,r, n);
}

// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float fOpUnionStairs(float a, float b, float r, float n) {
    float s = r/n;
    float u = b-r;
    return min(min(a,b), 0.5 * (u + a + abs ((mod (u - a + s, 2. * s)) - s)));
}

// We can just call Union since stairs are symmetric.
float fOpIntersectionStairs(float a, float b, float r, float n) {
    return -fOpUnionStairs(-a, -b, r, n);
}

float fOpDifferenceStairs(float a, float b, float r, float n) {
    return -fOpUnionStairs(-a, b, r, n);
}


// Similar to fOpUnionRound, but more lipschitz-y at acute angles
// (and less so at 90 degrees). Useful when fudging around too much
// by MediaMolecule, from Alex Evans' siggraph slides
float fOpUnionSoft(float a, float b, float r) {
    float e = max(r - abs(a - b), 0.);
    return min(a, b) - e*e*0.25/r;
}


// produces a cylindical pipe that runs along the intersection.
// No objects remain, only the pipe. This is not a boolean operator.
float fOpPipe(float a, float b, float r) {
    return length(vec2(a, b)) - r;
}

// first object gets a v-shaped engraving where it intersect the second
float fOpEngrave(float a, float b, float r) {
    return max(a, (a + r - abs(b))*sqrt(0.5));
}

// first object gets a capenter-style groove cut out
float fOpGroove(float a, float b, float ra, float rb) {
    return max(a, min(a + ra, rb - abs(b)));
}

// first object gets a capenter-style tongue attached
float fOpTongue(float a, float b, float ra, float rb) {
    return min(a, max(a - ra, abs(b) - rb));
}



// ---------------------------------------------------------------------- //

float bumpAmt = 0.0;
float strobeAmt = 10.0;
float globalTime = 0.0;
float darkAmt = 0.0;
float lightAmt = 0.0;
vec4 _m;

float noise1(vec2 uv){
    return texture2D(iChannel1, uv).x;   
}
float noise1s(vec2 uv){
    return noise1(uv)*2.0 - 1.0;
}
vec4 getMusic(){
    return _m;
}
vec4 blinnPhong(vec3 l, vec3 p, vec3 n, vec3 c)
{
    // Diffuse
    vec4 diff = vec4(c * max(dot(n, (l - p) / length(l - p)), 0.0) * 1., 1.0);
    // Specular
    return diff + vec4(c * 
                   pow(max(
                       dot(n, normalize(((l - p) / length(l - p)) 
                                        + normalize(-p))), 0.0), 64.0) * 1., 1.0);
}


float object(vec3 p, float r)
{
    //pMirrorOctant(p.xz, vec2(.55));
    //return fOctahedron(p, r, 22.);
    //return fTruncatedIcosahedron(p, r, 52.0);

    return fOpUnionRound(
        fDodecahedron(p, r, 32.0),
        fTorus(p, .1, .5),
        .1);
}


vec3 calcNormal(vec3 P, float r){
    vec3 eps = vec3(0.001, 0.0, 0.0);

    #define map(v) vec2(object(v, r))

    vec3 N = vec3(map(P+eps.xyy).x - map(P-eps.xyy).x,
                  map(P+eps.yxy).x - map(P-eps.yxy).x,
                  map(P+eps.yyx).x - map(P-eps.yyx).x);
    #undef map
    return normalize(N);
}

vec2 intersect(vec3 ro, vec3 rd, 
               out vec3 pOut, out vec3 nOut)
{
    vec3 xf = vec3(0.0, 0.0, 0.);
    float r = mix(.4,
                  getMusic().x * getMusic().x*.9,
                  .5*bumpAmt);
    
    // Matrix constructed as inverse.
    mat4 rot = rotationMatrix(vec3(0,0,1), 3.1415/2.);
    // Translation after roation
    rot[3].xyz = vec3(-1., 1., 0.);
            
    // ray parameter
    float t = 0.0;
    
    r*=0.75;
    
    for (float i = 0.0; i < 40.0; i += 1.) {
        vec3 p = ro + rd * t;
        vec4 pCyl = rot*vec4(p,1.0);
        float dd = object(p-xf, r);//sdSphere(p-xf, r);
        
        if (abs(dd) < 0.01) { 
            pOut = p;
            nOut = normalize(p-xf);
            nOut = calcNormal(p, r);
            return vec2(t, 1.0);
        }
        
        float dc = sdCylinder((pCyl).xyz, vec3(0.,0.,.1));
        if (abs(dc) < 0.01) {
            pOut = (vec4(p,1.)*rot).xyz;
            nOut = vec3(0);//normalize(p-xf);
            return vec2(t, 2.0);
        }
        
        float dc2 = sdCylinder(vec3(0.,0.,.4+ getMusic().z)+(pCyl).xyz, vec3(0.,0.,.1));
        if (abs(dc2) < 0.01) {
            pOut = (vec4(p,1.)*rot).xyz;
            nOut = vec3(0); //normalize(p-xf);
            return vec2(t, 2.0);
        }
        
        float dc3 = sdCylinder(vec3(0.,0.,-.4- getMusic().z)+(pCyl).xyz, vec3(0.,0.,.1));
        if (abs(dc3) < 0.01) {
            pOut = (vec4(p,1.)*rot).xyz;
            nOut = vec3(0); //normalize(p-xf);
            return vec2(t, 2.0);
        }

        t += min(min(min(abs(dd), abs(dc))
                 ,dc2)
                 ,dc3)-.008;
    }
    
    nOut = vec3(0);
    pOut = vec3(0);
    
    return vec2(0.0, 0.0);
}
// Specialized intersector just for tubes to make secondary rays
// super fast.
vec2 intersectTubes(vec3 ro, vec3 rd)
{
    mat4 rot = rotationMatrix(vec3(0,0,1), 3.1415/2.);
    rot[3].xyz = vec3(-1., 1., 0.);
    float t = 0.0;
    
    vec3 offsetPlus = vec3(0.,0.,.4+ getMusic().z);
    vec3 offsetMinus = vec3(0.,0.,-.4- getMusic().z);
    
    for (float i = 0.0; i < 10.0; i += 1.) {
        vec3 p = ro + rd * t;
        vec3 pCyl = (rot*vec4(p,1.0)).xyz;
        float dc = sdCylinder(pCyl, vec3(0.,0.,.1));        
        dc = min(dc, sdCylinder(offsetPlus+pCyl, vec3(0.,0.,.1)));
        dc = min(dc, sdCylinder(offsetMinus+pCyl, vec3(0.,0.,.1)));
        if (abs(dc) < 0.01) {
            return vec2(t, 2.0);
        }

        t += dc -.001;
    }
        
    return vec2(0.0, 0.0);
}

vec4 getColor2(vec3 ro, vec3 rd)
{
    // For speed, only accept hits from the tubes.
    return mix(vec4(0),
               vec4(1),
               vec4(intersectTubes(ro, rd).y == 2.0));
}
vec4 getColor(vec3 ro, vec3 rd, vec3 offset) {
    // Intersection
    vec3 p, n;
    vec2 hit;
    
    // hit.x = ray parameter
    // hit.y = object id
    hit = intersect(ro, rd, p, n);
    if (hit.x == 0.0) {
        return vec4(.5*textureCube(iChannel0, rd-offset)
                   +.5*textureCube(iChannel0, rd+offset));
    }
    vec4 normColor = vec4(.5 * n + .5, 1.0);
    vec4 hitColor = vec4(.0,.0,.0, 1.0);
    

    #if 0
    vec3 pp = vec3(0.0);
    vec3 nn = vec3(0.0);
    vec3 eps = vec3(0.001, 0.0, 0.0);
    #define map(RD) intersect(ro, RD, pp, nn)
    
    map(rd+eps.xyy);
    p = p + pp;
    n = n + nn;
    map(rd-eps.xyy);
    p = p + pp;
    n = n + nn;
    
    map(rd+eps.yxy);
    p = p + pp;
    n = n + nn;
    map(rd-eps.xyy);
    p = p + pp;
    n = n + nn;
    
    map(rd+eps.yxy);
    p = p + pp;
    n = n + nn;
    map(rd-eps.yxy);
    p = p + pp;
    n = n + nn;
    
    p /= 7.;
    n /= 7.;
    
    #undef map
    #endif
    
    
    
    // Lighting
    
    //return hitColor;
    hitColor += .75*vec4(textureCube(iChannel0, reflect(rd, n)).rgb, 1.0);
    
    hitColor += .8*getColor2(p + n*.01, reflect(rd,n));
    
    float flake = texture2D(iChannel1, n.yz*3.0).x;
    flake += .1*texture2D(iChannel1, n.xz).x;
    flake = 0.;
    float fresnelDarken = dot(-rd,n);
    vec3 c = .5*vec3(1.,1,1.);
    fresnelDarken *= fresnelDarken;
    hitColor += fresnelDarken*flake*blinnPhong(vec3(3,1,0), p, n, c);
    hitColor += fresnelDarken*flake*blinnPhong(vec3(-3,1,0), p, n, c);
    hitColor += fresnelDarken*flake*blinnPhong(vec3(0,1,0), p, n, c);
    
    return mix(hitColor,
               vec4(1),
               vec4(hit.y != 1.0));    
}
vec4 strobe(vec4 color) 
{
    return 1.7*color 
         * 1.5 * vec4(getMusic().x*.6, getMusic().x*.4, .3*getMusic().x, 1.0);
}
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    hg_sdf_init();
    // XXX: should be vec2(0.), but requres vec2(1.,0.)
    _m = texture2D(iChannel2, vec2(1.0,0.));
    // Mouse controls
    vec2 mo = iMouse.xy/iResolution.xy;
    mo -= vec2(.5,.5);

    globalTime = iChannelTime[2] + .4;
    // drop
    float fadeIntro = mix(0.0, 1.0, globalTime/5.85);
    float fadeIn = smoothstep(60.0, 61.0, globalTime);
    float fadeSpin = smoothstep(17.0, 20.5, globalTime);
    float fadeBump = smoothstep(20.0, 30.5, globalTime);
    float fadeCamSpaz = fadeBump;
    float fadeFirstBreak = smoothstep(28.0, 30.5, globalTime);
    bumpAmt = .1+fadeBump;
    
    fadeSpin = fadeSpin*(.9+getMusic().x/45.);
        
    // Fade in: 6s
    if (fadeIntro < 1.0) {
        fragColor = vec4(vec3(fadeIntro-.15), 1.0);
        return;
    }
    
    if (globalTime < 16.25) {
        bumpAmt = 1.0;   
    }

    if (globalTime > 28.0 && globalTime < 39.2) {
        fadeCamSpaz = 0.0;
        fadeSpin = 0.0;
        bumpAmt = mix(0.2,0.6, (globalTime-28.0)/(39.2-28.0));
    }
        
    if (globalTime > 228.) {
        fragColor = vec4(vec3(1.-(globalTime-228.)/(233.-228.)),1.0);
    }
    
    // Calm, looing up from below
    if (globalTime > 27.0 && globalTime < 59.5) {
        float s = globalTime/59.5;
        fadeCamSpaz = mix(0.0, fadeCamSpaz, s);
        fadeSpin  = mix(0.0, fadeSpin, s);
        bumpAmt = mix(0.0, bumpAmt, s);
        
        // Build up, looking down
        if (globalTime < 39.2) {
            mo += vec2(.5, 1.0);
        }
        if (globalTime > 59.5 && globalTime < 60.) {
            mo += vec2(.5, 1.0);
        }
    }
    
    // First Drop, looking up from below
    if (globalTime > 61.0 && globalTime < 85.) {
        mo += vec2(.5, 1.0);
    }

    // 1:25 - grr, looking up from below
    if (globalTime > 126.8 && globalTime < 128.) {
        mo += vec2(.5, 1.0)+.4*getMusic().x;
        //darkAmt = 1.0;
    }
    
    // faster at 116
    // 138 - drop out
    // 148 - build
    // 160 - drop in
    if (( globalTime > 0.5 && fadeIntro >= 1.0 && globalTime < 6.1)
       ||(globalTime > 16.25 && fadeIn < 1.0) 
       ||(globalTime > 104.5 && globalTime < 160.5)
    ) {
        strobeAmt = fadeIn = getMusic().x;
        float s = (globalTime-104.5)/(160.0-104.5);
        fadeCamSpaz = mix(0.0, fadeCamSpaz, s);
        //fadeSpin  = 0.5; //mix(0., 0.1, globalTime/110.); //mix(0.0, fadeSpin, s);
        fadeSpin  = mix(0.0, fadeSpin, s);
        bumpAmt = mix(0.0, bumpAmt, s);

        if ((globalTime > 60. && globalTime < 116.) 
            || (globalTime > 138. && globalTime < 140.5)
            || (globalTime > 143.5 && globalTime < 146.5 )
        ) {
            mo += vec2(.5, 1.0);
        }
    }
    
    if (globalTime > 160.5 && globalTime < 204.) {
        darkAmt = mix(0.3, 1.0, (globalTime-160.5)/(204.-160.5));
        if (globalTime < 170. || globalTime > 171.5)
            mo += vec2(.5, 1.0);
    }
    
    if ((globalTime > 203. && globalTime < 320.))
    {
        strobeAmt = fadeIn = getMusic().x;
        float s = (globalTime-204.)/(220.0-204.);
        fadeCamSpaz = mix(0.0, fadeCamSpaz, s);
        fadeSpin  = getMusic().x/45.; //mix(0.0, fadeSpin, s);
        bumpAmt = mix(0.0, bumpAmt, s);
        lightAmt = mix(0.0, 1.0, s);
        darkAmt = mix(1.0, 0.0, s);
    }
    
    if (( globalTime > 11.0 && globalTime < 11.2))
    {
        strobeAmt = fadeIn = getMusic().x;
    }
        
    
    // Camera, spherical coords
    float an1 = -6.2831*mo.x + 1.55 
                + mix(1.0,
                      (.1*globalTime*10.0 + noise1(vec2(globalTime)*.2)),
                      fadeSpin)
            ;
    float an2 = clamp(-.0  + 1.5*mo.y, 0.3, 3.35)
                + mix(mix(0.0, noise1(vec2(globalTime)*.2), fadeCamSpaz),
                      0.0, fadeFirstBreak);
    vec3 ro = normalize(vec3(sin(an2)*cos(an1), cos(an2)-0.5, sin(an2)*sin(an1)));
    
    vec2 uv = (2.*(fragCoord.xy+.5) - iResolution.xy) / iResolution.y;
    
    // Camera Transform
    vec3 ww = normalize(vec3(0.,0.,0.) - ro);
    vec3 uu = normalize(cross(vec3(0.0,1.0,0.0), ww));
    vec3 vv = normalize(cross(ww,uu));
    // Ray Direction
    vec3 rd = normalize(uv.x*uu + uv.y*vv + 1.4*ww);
    
    vec3 offset;
    {
        float m = getMusic().x*.02;
        float n1 = noise1s(vec2(globalTime));
        float n2 = noise1s(vec2(globalTime*7.));
        float n3 = noise1s(vec2(globalTime*3.));
        offset = vec3(m*n1,m*n2,m*n3);   
    }
    
    vec4 c = getColor(ro, rd, offset);
    
    #if 1
    // --------------------------------------------------------------- //
    // Brute force AA
    // --------------------------------------------------------------- //
    vec2 epSize = .5 / iResolution.xy;
    //epSize += epSize*.5 * noise1(epSize);
    vec3 eps = vec3(epSize.y, 0.0, -epSize.y);
    vec3 eps2 = eps*.75;
    #define map(RD) getColor(ro, RD, offset)
    c = (c
      + map(rd+eps.xyy) + map(rd-eps.xyy)
      + map(rd+eps.yxy) + map(rd-eps.yxy)
      + map(rd+eps.xxy) + map(rd-eps.xxy)
      + map(rd+eps.xzy) + map(rd-eps.xzy)
      ) / (8. + 1.);
    #undef map
    // --------------------------------------------------------------- //
    #endif
    
    c.r = mix(c.r,
              getColor(ro, rd 
                           +vec3(noise1(uv+sin(globalTime)*1e7)*.04 
                                      +(getMusic().x*.12-.066)),
                       offset).r,
              float(  (globalTime > 126.8 && globalTime < 128.)
                    || (globalTime > 80.5 && globalTime < 83.)
                    || (globalTime > 92.0 && globalTime < 93.5)
                    || (globalTime > 103.5 && globalTime < 105.)
                    || (globalTime > 138.0 && globalTime < 139.)
                    || (globalTime > 158. && globalTime < 160.)
                    || (globalTime > 180.5 && globalTime < 182.5)
                    || (globalTime > 191.5 && globalTime < 193.0)
                    || (globalTime > 203. && globalTime < 204.0)
                    || (globalTime > 59.5 && globalTime < 61.) 
                    )
             );

    
    c.b += mix(.0,  getMusic().x*.5, fadeIn);
    c.rgb = pow(c.rgb, vec3(1.3));
    fragColor = mix(strobe(c), c, 
                    1.0 - 
                    fadeIn*strobeAmt * 
                    (.5+.5*sin(2.*globalTime)));
    fragColor = mix(fragColor,
                    fragColor - .7*(length(uv)*length(uv)),
                    darkAmt);
    fragColor = mix(fragColor,
                    fragColor + .7*(length(uv)*length(uv)),
                    lightAmt);    
    fragColor = fragColor - .02*(length(uv)*length(uv)*length(uv)*length(uv));
}



