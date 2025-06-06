/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef ParticleInstancingRenderer_h
#define ParticleInstancingRenderer_h

#include <vector>
#include <DebugStream.h>

#ifndef VTK_IMPLEMENT_MESA_CXX
  #if defined(__APPLE__) && (defined(VTK_USE_CARBON) || defined(VTK_USE_COCOA))
    #include <OpenGL/gl.h>
  #else
    #if defined(_WIN32)
       #include <windows.h>
    #endif
    #include <GL/glew.h>
  #endif
#endif



class ParticleInstancingRenderer {
  public:
  ParticleInstancingRenderer();
  ~ParticleInstancingRenderer();

    void Initialize();
    bool IsSupported();

    void SetQualityLevel(int level);

    void AddParticle(const double* xyz, const double radius, const unsigned char* rgb);
    void ClearParticles();
    size_t NumParticles() const;

    void Render();
  
  private:
  
    void BuildShaders();
    void BuildSphereGeometryVBOs();

    void GenerateAndBuildTBO();
    void BuildTBO(const GLuint tbo, const GLuint tex, size_t tbo_size, GLenum usage,  GLenum internal_format);

  GLboolean CopyParticleDataToGpuBuffers(size_t start, size_t count,
                                         GLuint tbo_position_radius, GLuint tex_position_radius,
                                         GLuint tbo_color, GLuint tex_color);


    bool is_initialized;
    bool extensions_supported;

    GLuint program_instancing;
    size_t instanced_batch_size;
    
    GLuint tbo_position_radius_batches[2];
    GLuint tex_position_radius_batches[2];

    GLuint tbo_color_batches[2];
    GLuint tex_color_batches[2];
    
    int quality_level;
    std::vector<class SphereGeometryVBO*> sphere_geometry_vbos;

    std::vector<float> particle_position_radius;
    std::vector<unsigned char> particle_color;
};

#endif
