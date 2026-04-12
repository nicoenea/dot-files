// Original by localthunk (https://www.playbalatro.com)
// Modified for better terminal text contrast
 
// Configuration (modify these values to change the effect)
#define SPIN_ROTATION -2.0
#define SPIN_SPEED 7.0
#define OFFSET vec2(0.0)
#define COLOUR_1 vec4(0.85, 0.2, 0.2, 1.0)       // Vivid red (matches Balatro)
#define COLOUR_2 vec4(0.0, 0.612, 1.0, 1.0)     // Vivid blue (156/255, matches Balatro)
#define COLOUR_3 vec4(0.0, 0.0, 0.0, 1.0)       // Pure black background (matches Balatro)
#define CONTRAST 2.5
#define LIGTHING 0.35                            // Increased lighting for more detail
#define SPIN_AMOUNT 0.25
#define PIXEL_FILTER 745.0
#define SPIN_EASE 1.0
#define PI 3.14159265359
#define IS_ROTATE false
 
vec4 effect(vec2 screenSize, vec2 screen_coords) {
    float pixel_size = length(screenSize.xy) / PIXEL_FILTER;
    vec2 uv = (floor(screen_coords.xy*(1./pixel_size))*pixel_size - 0.5*screenSize.xy)/length(screenSize.xy) - OFFSET;
    float uv_len = length(uv);
 
    float speed = (SPIN_ROTATION*SPIN_EASE*0.2);
    if(IS_ROTATE){
       speed = iTime * speed;
    }
    speed += 302.2;
    float new_pixel_angle = atan(uv.y, uv.x) + speed - SPIN_EASE*20.*(1.*SPIN_AMOUNT*uv_len + (1. - 1.*SPIN_AMOUNT));
    vec2 mid = (screenSize.xy/length(screenSize.xy))/2.;
    uv = (vec2((uv_len * cos(new_pixel_angle) + mid.x), (uv_len * sin(new_pixel_angle) + mid.y)) - mid);
 
    uv *= 30.;
    speed = iTime*(SPIN_SPEED);
    vec2 uv2 = vec2(uv.x+uv.y);
 
    for(int i=0; i < 5; i++) {
        uv2 += sin(max(uv.x, uv.y)) + uv;
        uv  += 0.5*vec2(cos(5.1123314 + 0.353*uv2.y + speed*0.131121),sin(uv2.x - 0.113*speed));
        uv  -= 1.0*cos(uv.x + uv.y) - 1.0*sin(uv.x*0.711 - uv.y);
    }
 
    float contrast_mod = (0.25*CONTRAST + 0.5*SPIN_AMOUNT + 1.2);
    float paint_res = min(2., max(0.,length(uv)*(0.035)*contrast_mod));
    float c1p = max(0.,1. - contrast_mod*abs(1.-paint_res));
    float c2p = max(0.,1. - contrast_mod*abs(paint_res));
    float c3p = 1. - min(1., c1p + c2p);
    float light = (LIGTHING - 0.2)*max(c1p*5. - 4., 0.) + LIGTHING*max(c2p*5. - 4., 0.);
 
    // Create the base effect with more visible colors
    vec4 baseEffect = (0.3/CONTRAST)*COLOUR_1 + (1. - 0.3/CONTRAST)*(COLOUR_1*c1p + COLOUR_2*c2p + vec4(c3p*COLOUR_3.rgb, c3p*COLOUR_1.a)) + light;
 
    // Light darkening to maintain detail visibility
    baseEffect.rgb *= 1; // Only reduce brightness by 15%
 
    return baseEffect;
}
 
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / iResolution.xy;
 
    vec4 terminalColor = texture(iChannel0, uv);
    vec4 shaderEffect = effect(iResolution.xy, uv * iResolution.xy);
 
    // Balanced blending - prioritize text readability but keep background vibrant
    float hasContent = step(0.01, terminalColor.a);
 
    // Less aggressive text blending to preserve background beauty
    float textContrast = hasContent * 0.75 + (1.0 - hasContent) * 0.0;
 
    fragColor = mix(shaderEffect, terminalColor, textContrast);
}