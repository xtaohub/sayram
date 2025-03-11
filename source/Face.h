/*
 * File:        Face.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * 
 * Date:        11/13/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef FACE_H_
#define FACE_H_

#include "common.h"

enum Direction {XPOS, XNEG, YPOS, YNEG, ZPOS, ZNEG}; 

typedef std::array<Point, 4> Face_Vertices;

class Face{
  public:
    static const int NT = 4; // number of tetrahedrons per face. 
    const Point& v(int i) const { return vs_[i]; } // vertices
    int v_size() const { return vs_.size(); }                                                   
    double area() const { return area_; }
    const Vector3& n() const { return n_; } // normal vector, depending on the face direction
    Direction dir() const { return dir_; }

    /*
     * 
     * vs of the face should be in clockwise order 
     * when viewed from the center of the cell.
     * This should be ensured by whatever caller of set_vs_dir.
     * 
     */
    void set_vs_dir(const Face_Vertices& vs, Direction dir) {
      vs_ = vs; 
      area_ = (vs_[1] - vs_[0]).norm() * (vs_[2] - vs_[1]).norm(); 
      dir_ = dir; 

      switch (dir) {
        case XPOS: // Face perpendicular to positive x-axis (dy * dz)
          n_ = {1, 0, 0};
          break;

        case XNEG:
          n_ = {-1, 0, 0};
          break; 

        case YPOS: // Face perpendicular to y-axis (dx * dz)
          n_ = {0, 1, 0};
          break;

        case YNEG:
          n_ = {0, -1, 0};
          break;

        case ZPOS: // Face perpendicular to z-axis (dx * dy)
          n_ = {0, 0, 1}; 
          break;

        case ZNEG:
          n_ = {0, 0, -1}; 
          break;

        default: // incorrect dir
          n_ = {0, 0, 0}; 
      }
    }

  private:

    Face_Vertices vs_; 
    double area_; 
    Vector3 n_; // the unit normal vector
    Direction dir_; 
};

#endif /* FACE_H */

