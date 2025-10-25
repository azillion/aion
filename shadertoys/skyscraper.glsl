// one skyscraper: doors + no first-floor windows + fog
// static camera; no time; single pass

#define MAX_STEPS 96
#define MAX_DIST  100.0
#define SURF_EPS  0.0018

// building
const vec3  BLDG_SIZE   = vec3(2.5, 7.5, 2.5);   // half-extents (x,y,z)
const float LOBBY_H     = 1.4;                   // first floor height (world units)

// (doors and lobby inset removed)
// ground floor (windowless band height)
const float GROUND_FLOOR_H = 1.35;               // windowless band height from ground

// camera
const vec3  CAM_POS     = vec3(8.0, 4.0, 10.0);
const vec3  CAM_TGT     = vec3(0.0, 3.0, 0.0);
const float FOV         = 1.1;
const float FOG_DENS    = 0.05;

// materials
const vec3  ALB_CONC    = vec3(0.18,0.19,0.21);
const vec3  ALB_GROUND  = vec3(0.08,0.085,0.09);
const vec3  ALB_GLASS   = vec3(0.05,0.07,0.09);
const vec3  LOBBY_COL   = vec3(1.0, 0.92, 0.78); // warm lobby
const float LOBBY_INT   = 1.6;

// window system (auto-fit) — Art Deco: narrow, tall windows with slim gaps
// world units (no cells)
const float WIN_MARGIN  = 0.02;
const float WIN_GAP_X   = 0.03;
const float WIN_GAP_Y   = 0.05;
const float WIN_MIN_X   = 0.22;
const float WIN_MIN_Y   = 0.50;
const float WIN_FRAME   = 0.12;
const float WIN_EDGE_PAD = 0.04;              // normalized edge guard on faces

// Art Deco tiered setbacks (fractions of height and plan scale)
const float TIER1_Y_FRAC   = 0.53;
const float TIER2_Y_FRAC   = 0.78;
const float TIER1_XZ_SCALE = 0.85;
const float TIER2_XZ_SCALE = 0.68;
const float LIT_RATIO   = 0.65;
const vec3  LIGHT_COL   = vec3(1.00, 0.88, 0.65);
const float LIGHT_MIN   = 0.7;
const float LIGHT_MAX   = 2.6;
const float REVEAL_DEPTH = 0.04;                 // window inset depth

float hash21(vec2 p){ return fract(sin(dot(p, vec2(127.1,311.7)))*43758.5453); }
mat2 rot(float a){ float s=sin(a), c=cos(a); return mat2(c,-s,s,c); }

float sdBox(vec3 p, vec3 b){
    vec3 q = abs(p) - b;
    return length(max(q,0.0)) + min(max(q.x, max(q.y,q.z)), 0.0);
}
float sdPlaneY(vec3 p){ return p.y; }
float sdRoundBox(vec3 p, vec3 b, float r){
    vec3 q = abs(p) - (b - vec3(r));
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}

// pick the active Art Deco tier for a given world-space y
void getDecoTier(float y, out vec3 b, out vec3 c, out int tierIdx){
    float yTop = BLDG_SIZE.y;
    float y0 = -yTop;
    float y1 = -yTop + 2.0*yTop*TIER1_Y_FRAC;
    float y2 = -yTop + 2.0*yTop*TIER2_Y_FRAC;
    float y3 = yTop;
    float scale; float bottom; float top; int idx;
    if(y < y1){ scale=1.0; bottom=y0; top=y1; idx=0; }
    else if(y < y2){ scale=TIER1_XZ_SCALE; bottom=y1; top=y2; idx=1; }
    else{ scale=TIER2_XZ_SCALE; bottom=y2; top=y3; idx=2; }
    b = vec3(BLDG_SIZE.x*scale, 0.5*(top-bottom), BLDG_SIZE.z*scale);
    c = vec3(0.0, 0.5*(top+bottom), 0.0);
    tierIdx = idx;
}

// ---------- face coords via sdf-contact axis ----------
struct Face2D { vec2 xy; vec2 ext; vec3 fid; };
Face2D boxFace2D(vec3 p, vec3 b){
    vec3 q = abs(p) - b;
    vec3 m = step(q.yzx, q.xyz) * step(q.zxy, q.xyz);
    Face2D F;
    if(m.x>0.5){ F.xy=vec2(p.z,p.y); F.ext=vec2(b.z,b.y); F.fid=vec3(sign(p.x),0,0); }
    else if(m.y>0.5){ F.xy=vec2(p.x,p.z); F.ext=vec2(b.x,b.z); F.fid=vec3(0,sign(p.y),0); }
    else            { F.xy=vec2(p.x,p.y); F.ext=vec2(b.x,b.y); F.fid=vec3(0,0,sign(p.z)); }
    return F;
}

// axis packer
void fitAxis(float E, float M, float G, float Wmin,
             out int N, out float W, out float g, out float S){
    float usable = max(0.0, 2.0*E - 2.0*M);
    int n = int(floor( (usable + G) / (Wmin + G) ));
    n = max(n, 0);
    if(n==0 && usable >= Wmin) n=1;
    if(n==0){ N=0; W=0.0; g=0.0; S=0.0; return; }
    if(n==1){ W=usable; g=0.0; }
    else{
        float Wcand = (usable - G*float(n-1)) / float(n);
        if(Wcand < Wmin){
            for(int k=0;k<16;k++){ n--; if(n<=1) break; Wcand=(usable - G*float(n-1))/float(n); if(Wcand>=Wmin) break; }
            n=max(n,1);
            if(n==1){ W=usable; g=0.0; } else { W=Wcand; g=G; }
        }else{ W=Wcand; g=G; }
    }
    float total = float(n)*W + float(max(n-1,0))*g;
    float left  = -E + M + 0.5*W + 0.5*(usable - total);
    N=n; S=left;
}

// returns paneMask (0/1 glass), lightVal (0..amp; 0 means OFF)
void windowAutoParams(vec3 p, vec3 b,
                      float margin, float gapX, float gapY,
                      float minWX, float minWY, float frame,
                      out float paneMask, out float lightVal)
{
    paneMask = 0.0; lightVal = 0.0;

    // tier-local coordinates and extents
    vec3 tb, tc; int tidx; getDecoTier(p.y, tb, tc, tidx);
    vec3 pl = p - tc;
    Face2D F = boxFace2D(pl, tb);
    if(abs(F.fid.y) > 0.5) return; // no roof/ground

    // no windows in the first floor band (only for base tier)
    float groundY = -F.ext.y;
    if(tidx == 0 && F.xy.y < groundY + GROUND_FLOOR_H) return;

    // edge guard: avoid windows too close to vertical edges
    float ex = (F.ext.x - abs(F.xy.x)) / F.ext.x;
    float ey = (F.ext.y - abs(F.xy.y)) / F.ext.y;
    float edgeNorm = min(ex, ey);
    if(edgeNorm < WIN_EDGE_PAD) return;

    int nx, ny; float wx, wy, gx, gy, sx, sy;
    fitAxis(F.ext.x, margin, gapX, minWX, nx, wx, gx, sx);
    fitAxis(F.ext.y, margin, gapY, minWY, ny, wy, gy, sy);
    if(nx<=0 || ny<=0) return;

    float kx = (F.xy.x - sx) / (wx + gx);
    float ky = (F.xy.y - sy) / (wy + gy);
    int ix = int(floor(kx + 1e-4));
    int iy = int(floor(ky + 1e-4));
    if(ix<0 || ix>=nx || iy<0 || iy>=ny) return;

    float cx = sx + float(ix)*(wx+gx);
    float cy = sy + float(iy)*(wy+gy);
    vec2  d  = abs(F.xy - vec2(cx, cy));
    vec2  halfWin = 0.5*vec2(wx, wy);

    // inside the window rect at all?
    // require a minimum headroom to avoid black seams at tier roofs
    if(any(greaterThan(d, halfWin))) return;
    float topY = F.ext.y;
    float headroom = topY - (cy + halfWin.y);
    // allow tighter fit on top tier to reach near the crown
    float minHead = (tidx==2) ? 0.005 : max(0.02, frame*wy);
    if(headroom < minHead) return;

    // glass pane after frame cut
    vec2 inner = halfWin * (1.0 - frame*2.0);
    float pane = step(d.x, inner.x) * step(d.y, inner.y);

    // per-window ON/OFF + intensity
    vec2 key = vec2(float(ix), float(iy)) + F.fid.xy*7.31;
    float r0 = hash21(key);
    float r1 = hash21(key + 13.17);
    float on  = step(r0, LIT_RATIO);
    float amp = mix(LIGHT_MIN, LIGHT_MAX, r1);

    paneMask = pane;
    lightVal = pane * on * amp;
}

// fast pane mask at arbitrary point (for SDF reveal)
float windowPaneMaskAt(vec3 p, vec3 b,
                       float margin, float gapX, float gapY,
                       float minWX, float minWY, float frame)
{
    vec3 tb, tc; int tidx; getDecoTier(p.y, tb, tc, tidx);
    vec3 pl = p - tc;
    Face2D F = boxFace2D(pl,tb);
    if(abs(F.fid.y) > 0.5) return 0.0; // skip roof/ground

    float groundY = -F.ext.y;
    if(tidx == 0 && F.xy.y < groundY + GROUND_FLOOR_H) return 0.0; // no first-floor windows only on base tier

    // edge guard: avoid windows too close to vertical edges
    float ex = (F.ext.x - abs(F.xy.x)) / F.ext.x;
    float ey = (F.ext.y - abs(F.xy.y)) / F.ext.y;
    float edgeNorm = min(ex, ey);
    if(edgeNorm < WIN_EDGE_PAD) return 0.0;

    int nx, ny; float wx, wy, gx, gy, sx, sy;
    fitAxis(F.ext.x, margin, gapX, minWX, nx, wx, gx, sx);
    fitAxis(F.ext.y, margin, gapY, minWY, ny, wy, gy, sy);
    if(nx<=0 || ny<=0) return 0.0;

    float kx = (F.xy.x - sx) / (wx + gx);
    float ky = (F.xy.y - sy) / (wy + gy);
    int ix = int(floor(kx + 1e-5));
    int iy = int(floor(ky + 1e-5));
    if(ix<0 || ix>=nx || iy<0 || iy>=ny) return 0.0;

    float cx = sx + float(ix)*(wx+gx);
    float cy = sy + float(iy)*(wy+gy);
    vec2 d = abs(F.xy - vec2(cx, cy));
    vec2 halfWin = 0.5*vec2(wx, wy);
    if(any(greaterThan(d, halfWin))) return 0.0;
    float topY = F.ext.y;
    float headroom = topY - (cy + halfWin.y);
    float minHead = (tidx==2) ? 0.005 : max(0.02, frame*wy);
    if(headroom < minHead) return 0.0;

    vec2 inner = halfWin * (1.0 - frame*2.0);
    float pane = step(d.x, inner.x) * step(d.y, inner.y);
    return pane;
}

// (door masks removed)

// ---------- scene sdf (stacked shrinking boxes) ----------
struct Hit { float d; int mat; }; // mat: 0 ground, 1 building
Hit mapScene(vec3 p){
    float dGround = sdPlaneY(p);

    // stacked tiers: choose active tier by y and evaluate a single round box centered on that tier
    vec3 tb, tc; int tidx; getDecoTier(p.y, tb, tc, tidx);
    float dBuild = sdRoundBox(p - tc, tb, 0.06);

    // window reveal: inset glass areas slightly to create sill depth (push inward)
    float paneHere = windowPaneMaskAt(p, BLDG_SIZE,
                                      WIN_MARGIN, WIN_GAP_X, WIN_GAP_Y,
                                      WIN_MIN_X, WIN_MIN_Y, WIN_FRAME);
    float revealSDF = (paneHere > 0.0) ? REVEAL_DEPTH : 0.0;
    dBuild = max(dBuild, dBuild + revealSDF);

    Hit h;
    if(dBuild < dGround){ h.d = dBuild; h.mat = 1; }
    else                { h.d = dGround; h.mat = 0; }
    return h;
}

vec3 calcNormal(vec3 p){
    vec2 e = vec2(0.001, -0.001);
    return normalize(
        e.xyy * mapScene(p + e.xyy).d +
        e.yyx * mapScene(p + e.yyx).d +
        e.yxy * mapScene(p + e.yxy).d +
        e.xxx * mapScene(p + e.xxx).d
    );
}
float softShadow(vec3 ro, vec3 rd, float k){
    float res=1.0, t=0.02;
    for(int i=0;i<32;i++){
        float h = mapScene(ro + rd*t).d;
        res = min(res, k*h/t);
        t += clamp(h, 0.02, 0.6);
        if(h<0.0006 || t>40.0) break;
    }
    return clamp(res,0.0,1.0);
}
float ao(vec3 p, vec3 n){
    float occ=0.0, sca=1.0;
    for(int i=0;i<5;i++){
        float hr=0.03+0.12*float(i*i);
        float d = mapScene(p + n*hr).d;
        occ += (hr - d)*sca;
        sca *= 0.75;
    }
    return clamp(1.0 - 1.25*occ, 0.0, 1.0);
}

// camera
struct Ray { vec3 ro; vec3 rd; };
Ray makeCamera(vec2 uv, vec3 pos, vec3 ta, float fov){
    vec3 f = normalize(ta - pos);
    vec3 r = normalize(cross(vec3(0,1,0), f));
    vec3 u = cross(f, r);
    vec3 rd = normalize(r*uv.x + u*uv.y + f*fov);
    Ray ray; ray.ro=pos; ray.rd=rd; return ray;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord){
    vec2 R = iResolution.xy;
    vec2 uv = (fragCoord - 0.5*R) / R.y;

    Ray ray = makeCamera(uv, CAM_POS, CAM_TGT, FOV);

    float t=0.0; int mat=-1; int i;
    for(i=0;i<MAX_STEPS;i++){
        Hit h = mapScene(ray.ro + ray.rd*t);
        if(h.d < SURF_EPS || t > MAX_DIST){ mat=h.mat; break; }
        t += h.d;
    }

    vec3 col=vec3(0.0);
    float fogAmt = 1.0 - exp(-FOG_DENS * t);

    if(t<MAX_DIST && i<MAX_STEPS){
        vec3 p = ray.ro + ray.rd*t;
        vec3 n = calcNormal(p);

        vec3 sunDir = normalize(vec3(0.6,1.0,0.35));
        float ndl = max(dot(n,sunDir),0.0);
        float sh  = softShadow(p + n*0.02, sunDir, 12.0);
        float occ = ao(p,n);

        // base albedo — Art Deco: mostly stone with glass windows
        vec3 albedo = (mat==1)? ALB_CONC : ALB_GROUND;

        // windows (skip first floor)
        float paneMask=0.0, lightVal=0.0;
        if(mat==1){
            windowAutoParams(p, BLDG_SIZE,
                            WIN_MARGIN, WIN_GAP_X, WIN_GAP_Y,
                            WIN_MIN_X, WIN_MIN_Y, WIN_FRAME,
                             paneMask, lightVal);
            // now paneMask represents glass tiles within glass wall; keep albedo glass
        }

        // (doors removed)

        // lighting terms
        vec3 diff = albedo * ndl * sh;
        vec3 amb  = albedo * (0.22 + 0.22*occ);

        // spec + extra on glass (fresnel)
        float NoV = max(dot(n, -ray.rd), 0.0);
        float fres = pow(1.0 - NoV, 5.0);
        float spec = pow(max(dot(n, normalize(sunDir - ray.rd)), 0.0), 64.0) * sh;
        // stone base with glass tiles: modest base spec, bigger on panes
        spec *= mix(0.35, 0.75, paneMask) * (1.0 + fres*0.5);
        // (no door spec boost)

        // emissive: windows only
        vec3 em = LIGHT_COL * lightVal;

        // (no lobby panel emission)

        // wet base near ground
        float wet = smoothstep(0.02, 0.0, p.y) * smoothstep(1.2, 0.3, length(p.xz));
        vec3 albedoWet = mix(albedo, albedo*0.6, wet);
        float specWet  = mix(1.0, 2.0, wet);

        vec3 lit = amb + (albedoWet / max(albedo, vec3(1e-3))) * diff + spec*specWet + em;
        col = lit;
    }

    // fog
    vec3 fogCol = vec3(0.02, 0.03, 0.05);
    col = mix(col, fogCol, clamp(fogAmt,0.0,1.0));

    // vignette + gamma
    float vig = smoothstep(1.35, 0.25, length(uv*vec2(1.2,1.0)));
    col *= vig;
    col = pow(max(col,0.0), vec3(0.9));

    fragColor = vec4(clamp(col,0.0,1.0),1.0);
}
