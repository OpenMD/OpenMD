#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hydrodynamics/Mesh.hpp"

namespace OpenMD {

#define MSMSLINESIZE 180

  class MSMSFormat {
  public:
    MSMSFormat(const std::string& msmsRoot);
    ~MSMSFormat();
    Mesh* ReadShape();
    std::string title_;

  private:
    Mesh* mesh_;
    std::istream* faceFile_;
    std::istream* vertFile_;
    std::vector<Vector3d> vertices_;
    std::vector<Vector3d> normals_;
  };

  MSMSFormat::MSMSFormat(const std::string& msmsRoot) {
    std::string prefix = getPrefix(msmsRoot);  // just in case the user
                                               // gave us one of the
                                               // MSMS .face or .vert
                                               // file names

    std::string faceFileName = prefix + ".face";
    std::string vertFileName = prefix + ".vert";

    faceFile_ = new std::ifstream(faceFileName.c_str(),
                                  ifstream::in | ifstream::binary);

    if (faceFile_->fail()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "MSMSFormat: Cannot open facet file: %s\n",
               faceFileName.c_str());
      painCave.isFatal = 1;
      simError();
    }

    vertFile_ = new std::ifstream(vertFileName.c_str(),
                                  ifstream::in | ifstream::binary);

    if (vertFile_->fail()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "MSMSFormat: Cannot open vertex file: %s\n",
               vertFileName.c_str());
      painCave.isFatal = 1;
      simError();
    }
    mesh_ = new Mesh();

    return;
  }

  MSMSFormat::~MSMSFormat() {
    delete faceFile_;
    delete vertFile_;
    return;
  }

  Mesh* MSMSFormat::ReadShape() {
    char inbuf[MSMSLINESIZE];
    std::string line;

    vertFile_->clear();
    vertFile_->seekg(0);

    while (vertFile_->getline(inbuf, MSMSLINESIZE)) {
      // ignore comment lines:
      if (inbuf[0] != '#') {
        line = inbuf;
        StringTokenizer tokenizer(line);
        // do we have at least 3 tokens on this vertex line?
        if (tokenizer.countTokens() >= 3) {
          Vector3d vertex;
          // int l0fa, atomid, l;
          vertex[0] = tokenizer.nextTokenAsDouble();
          vertex[1] = tokenizer.nextTokenAsDouble();
          vertex[2] = tokenizer.nextTokenAsDouble();
          vertices_.push_back(vertex);
        }
        if (tokenizer.countTokens() >= 6) {
          Vector3d normal;
          normal[0] = tokenizer.nextTokenAsDouble();
          normal[1] = tokenizer.nextTokenAsDouble();
          normal[2] = tokenizer.nextTokenAsDouble();
          normals_.push_back(normal);                    
        }
        // l0fa = tokenizer.nextTokenAsInt();
        // atomid = tokenizer.nextTokenAsInt();
        // l = tokenizer.nextTokenAsInt();               
      }
    }

    faceFile_->clear();
    faceFile_->seekg(0);
    while (faceFile_->getline(inbuf, MSMSLINESIZE)) {
      // ignore comment lines:
      if (inbuf[0] != '#') {
        line = inbuf;
        StringTokenizer tokenizer(line);
        // do we have at least 3 tokens on this vertex line?
        if (tokenizer.countTokens() >= 3) {
          int v0, v1, v2;
          // int surftype, ana;
          v0 = tokenizer.nextTokenAsInt();
          v1 = tokenizer.nextTokenAsInt();
          v2 = tokenizer.nextTokenAsInt();
          // surftype = tokenizer.nextTokenAsInt();
          // ana = tokenizer.nextTokenAsInt();
          v0--;  // convert from 1-based indexing to 0-based indexing
          v1--;
          v2--;
          mesh_->add(vertices_[v0], vertices_[v1], vertices_[v2]);
        }
      }
    }
    std::cerr << "MSMS files have " << vertices_.size() << " vertices and ";
    std::cerr << mesh_->size() << " facets\n";
    return mesh_;
  }
}  // namespace OpenMD
