void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / iResolution.xy;
    vec4 original = texture(iChannel0, uv);
    
    // Add a red tint - this should be very obvious
    fragColor = vec4(original.rgb + vec3(0.2, 0.0, 0.0), original.a);
}