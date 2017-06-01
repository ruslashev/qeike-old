#include "utils.hh"

namespace shaders {

const char *light_vert = _glsl(
  attribute vec3 vertex_pos;
  attribute vec3 vertex_normal;
  uniform mat4 model;
  uniform mat4 view;
  uniform mat4 projection;
  varying vec3 frag_normal;
  varying vec3 frag_pos;
  void main() {
    frag_pos = vec3(model * vec4(vertex_pos, 1.0f));
    frag_normal = normalize(mat3(transpose(inverse(model))) * vertex_normal);
    gl_Position = projection * view * model * vec4(vertex_pos, 1.0);
  }
);

const char *light_frag = _glsl(
  uniform vec2 iResolution;
  uniform float iGlobalTime;
  uniform vec3 light_pos;
  uniform vec3 view_pos;
  uniform vec3 object_color;
  varying vec3 frag_normal;
  varying vec3 frag_pos;
  void main() {
    vec3 light_color = vec3(1, 1, 1);
    float ambient_strength = 0.1f;
    vec3 ambient = ambient_strength * light_color;

    vec3 light_dir = normalize(light_pos - frag_pos);
    float diff = max(dot(frag_normal, light_dir), 0.0);
    vec3 diffuse = diff * light_color;

    float specular_strength = 0.5f;
    vec3 view_dir = normalize(view_pos - frag_pos);
    vec3 reflect_dir = reflect(-light_dir, frag_normal);
    float spec = pow(max(dot(view_dir, reflect_dir), 0.0), 32);
    vec3 specular = specular_strength * spec * light_color;

    vec3 result = (ambient + diffuse + specular) * object_color;
    gl_FragColor = vec4(result, 1.0);
  }
);

const char *map_vert = _glsl(
  in vec3 Position;
  in vec2 TextureCoord;
  in vec2 LightmapCoord;
  uniform mat4 model;
  uniform mat4 view;
  uniform mat4 projection;
  out vec2 TexCoord0;
  out vec2 TexCoord1;
  void main()
  {
    gl_Position = projection * view * model * vec4(Position, 1.0);
    TexCoord0 = TextureCoord;
    TexCoord1 = LightmapCoord;
  }
);

const char *map_frag = _glsl(
  in vec2 TexCoord0;
  in vec2 TexCoord1;
  out vec4 FragColor;
  uniform sampler2D textureSampler;
  uniform sampler2D lightmapSampler;
  void main()
  {
    FragColor = texture(lightmapSampler, TexCoord1);//lightmap
  }
);

};

