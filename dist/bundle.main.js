/******/ (() => {
  // webpackBootstrap
  /******/ "use strict";
  /******/ var __webpack_modules__ = {
    /***/ "./node_modules/cannon-es/dist/cannon-es.js":
      /*!**************************************************!*\
  !*** ./node_modules/cannon-es/dist/cannon-es.js ***!
  \**************************************************/
      /***/ (
        __unused_webpack_module,
        __webpack_exports__,
        __webpack_require__
      ) => {
        __webpack_require__.r(__webpack_exports__);
        /* harmony export */ __webpack_require__.d(__webpack_exports__, {
          /* harmony export */ AABB: () => /* binding */ AABB,
          /* harmony export */ ArrayCollisionMatrix: () =>
            /* binding */ ArrayCollisionMatrix,
          /* harmony export */ BODY_SLEEP_STATES: () =>
            /* binding */ BODY_SLEEP_STATES,
          /* harmony export */ BODY_TYPES: () => /* binding */ BODY_TYPES,
          /* harmony export */ Body: () => /* binding */ Body,
          /* harmony export */ Box: () => /* binding */ Box,
          /* harmony export */ Broadphase: () => /* binding */ Broadphase,
          /* harmony export */ COLLISION_TYPES: () =>
            /* binding */ COLLISION_TYPES,
          /* harmony export */ ConeTwistConstraint: () =>
            /* binding */ ConeTwistConstraint,
          /* harmony export */ Constraint: () => /* binding */ Constraint,
          /* harmony export */ ContactEquation: () =>
            /* binding */ ContactEquation,
          /* harmony export */ ContactMaterial: () =>
            /* binding */ ContactMaterial,
          /* harmony export */ ConvexPolyhedron: () =>
            /* binding */ ConvexPolyhedron,
          /* harmony export */ Cylinder: () => /* binding */ Cylinder,
          /* harmony export */ DistanceConstraint: () =>
            /* binding */ DistanceConstraint,
          /* harmony export */ Equation: () => /* binding */ Equation,
          /* harmony export */ EventTarget: () => /* binding */ EventTarget,
          /* harmony export */ FrictionEquation: () =>
            /* binding */ FrictionEquation,
          /* harmony export */ GSSolver: () => /* binding */ GSSolver,
          /* harmony export */ GridBroadphase: () =>
            /* binding */ GridBroadphase,
          /* harmony export */ Heightfield: () => /* binding */ Heightfield,
          /* harmony export */ HingeConstraint: () =>
            /* binding */ HingeConstraint,
          /* harmony export */ JacobianElement: () =>
            /* binding */ JacobianElement,
          /* harmony export */ LockConstraint: () =>
            /* binding */ LockConstraint,
          /* harmony export */ Mat3: () => /* binding */ Mat3,
          /* harmony export */ Material: () => /* binding */ Material,
          /* harmony export */ NaiveBroadphase: () =>
            /* binding */ NaiveBroadphase,
          /* harmony export */ Narrowphase: () => /* binding */ Narrowphase,
          /* harmony export */ ObjectCollisionMatrix: () =>
            /* binding */ ObjectCollisionMatrix,
          /* harmony export */ Particle: () => /* binding */ Particle,
          /* harmony export */ Plane: () => /* binding */ Plane,
          /* harmony export */ PointToPointConstraint: () =>
            /* binding */ PointToPointConstraint,
          /* harmony export */ Pool: () => /* binding */ Pool,
          /* harmony export */ Quaternion: () => /* binding */ Quaternion,
          /* harmony export */ RAY_MODES: () => /* binding */ RAY_MODES,
          /* harmony export */ Ray: () => /* binding */ Ray,
          /* harmony export */ RaycastResult: () => /* binding */ RaycastResult,
          /* harmony export */ RaycastVehicle: () =>
            /* binding */ RaycastVehicle,
          /* harmony export */ RigidVehicle: () => /* binding */ RigidVehicle,
          /* harmony export */ RotationalEquation: () =>
            /* binding */ RotationalEquation,
          /* harmony export */ RotationalMotorEquation: () =>
            /* binding */ RotationalMotorEquation,
          /* harmony export */ SAPBroadphase: () => /* binding */ SAPBroadphase,
          /* harmony export */ SHAPE_TYPES: () => /* binding */ SHAPE_TYPES,
          /* harmony export */ SPHSystem: () => /* binding */ SPHSystem,
          /* harmony export */ Shape: () => /* binding */ Shape,
          /* harmony export */ Solver: () => /* binding */ Solver,
          /* harmony export */ Sphere: () => /* binding */ Sphere,
          /* harmony export */ SplitSolver: () => /* binding */ SplitSolver,
          /* harmony export */ Spring: () => /* binding */ Spring,
          /* harmony export */ Transform: () => /* binding */ Transform,
          /* harmony export */ Trimesh: () => /* binding */ Trimesh,
          /* harmony export */ Vec3: () => /* binding */ Vec3,
          /* harmony export */ Vec3Pool: () => /* binding */ Vec3Pool,
          /* harmony export */ WheelInfo: () => /* binding */ WheelInfo,
          /* harmony export */ World: () => /* binding */ World,
          /* harmony export */
        });
        /**
         * Records what objects are colliding with each other
         */
        class ObjectCollisionMatrix {
          /**
           * The matrix storage.
           */

          /**
           * @todo Remove useless constructor
           */
          constructor() {
            this.matrix = void 0;
            this.matrix = {};
          }
          /**
           * get
           */

          get(bi, bj) {
            let { id: i } = bi;
            let { id: j } = bj;

            if (j > i) {
              const temp = j;
              j = i;
              i = temp;
            }

            return i + "-" + j in this.matrix;
          }
          /**
           * set
           */

          set(bi, bj, value) {
            let { id: i } = bi;
            let { id: j } = bj;

            if (j > i) {
              const temp = j;
              j = i;
              i = temp;
            }

            if (value) {
              this.matrix[i + "-" + j] = true;
            } else {
              delete this.matrix[i + "-" + j];
            }
          }
          /**
           * Empty the matrix
           */

          reset() {
            this.matrix = {};
          }
          /**
           * Set max number of objects
           */

          setNumObjects(n) {}
        }

        /**
         * A 3x3 matrix.
         * Authored by {@link http://github.com/schteppe/ schteppe}
         */
        class Mat3 {
          /**
           * A vector of length 9, containing all matrix elements.
           */

          /**
           * @param elements A vector of length 9, containing all matrix elements.
           */
          constructor(elements = [0, 0, 0, 0, 0, 0, 0, 0, 0]) {
            this.elements = void 0;
            this.elements = elements;
          }
          /**
           * Sets the matrix to identity
           * @todo Should perhaps be renamed to `setIdentity()` to be more clear.
           * @todo Create another function that immediately creates an identity matrix eg. `eye()`
           */

          identity() {
            const e = this.elements;
            e[0] = 1;
            e[1] = 0;
            e[2] = 0;
            e[3] = 0;
            e[4] = 1;
            e[5] = 0;
            e[6] = 0;
            e[7] = 0;
            e[8] = 1;
          }
          /**
           * Set all elements to zero
           */

          setZero() {
            const e = this.elements;
            e[0] = 0;
            e[1] = 0;
            e[2] = 0;
            e[3] = 0;
            e[4] = 0;
            e[5] = 0;
            e[6] = 0;
            e[7] = 0;
            e[8] = 0;
          }
          /**
           * Sets the matrix diagonal elements from a Vec3
           */

          setTrace(vector) {
            const e = this.elements;
            e[0] = vector.x;
            e[4] = vector.y;
            e[8] = vector.z;
          }
          /**
           * Gets the matrix diagonal elements
           */

          getTrace(target = new Vec3()) {
            const e = this.elements;
            target.x = e[0];
            target.y = e[4];
            target.z = e[8];
            return target;
          }
          /**
           * Matrix-Vector multiplication
           * @param v The vector to multiply with
           * @param target Optional, target to save the result in.
           */

          vmult(v, target = new Vec3()) {
            const e = this.elements;
            const x = v.x;
            const y = v.y;
            const z = v.z;
            target.x = e[0] * x + e[1] * y + e[2] * z;
            target.y = e[3] * x + e[4] * y + e[5] * z;
            target.z = e[6] * x + e[7] * y + e[8] * z;
            return target;
          }
          /**
           * Matrix-scalar multiplication
           */

          smult(s) {
            for (let i = 0; i < this.elements.length; i++) {
              this.elements[i] *= s;
            }
          }
          /**
           * Matrix multiplication
           * @param matrix Matrix to multiply with from left side.
           */

          mmult(matrix, target = new Mat3()) {
            const A = this.elements;
            const B = matrix.elements;
            const T = target.elements;
            const a11 = A[0],
              a12 = A[1],
              a13 = A[2],
              a21 = A[3],
              a22 = A[4],
              a23 = A[5],
              a31 = A[6],
              a32 = A[7],
              a33 = A[8];
            const b11 = B[0],
              b12 = B[1],
              b13 = B[2],
              b21 = B[3],
              b22 = B[4],
              b23 = B[5],
              b31 = B[6],
              b32 = B[7],
              b33 = B[8];
            T[0] = a11 * b11 + a12 * b21 + a13 * b31;
            T[1] = a11 * b12 + a12 * b22 + a13 * b32;
            T[2] = a11 * b13 + a12 * b23 + a13 * b33;
            T[3] = a21 * b11 + a22 * b21 + a23 * b31;
            T[4] = a21 * b12 + a22 * b22 + a23 * b32;
            T[5] = a21 * b13 + a22 * b23 + a23 * b33;
            T[6] = a31 * b11 + a32 * b21 + a33 * b31;
            T[7] = a31 * b12 + a32 * b22 + a33 * b32;
            T[8] = a31 * b13 + a32 * b23 + a33 * b33;
            return target;
          }
          /**
           * Scale each column of the matrix
           */

          scale(vector, target = new Mat3()) {
            const e = this.elements;
            const t = target.elements;

            for (let i = 0; i !== 3; i++) {
              t[3 * i + 0] = vector.x * e[3 * i + 0];
              t[3 * i + 1] = vector.y * e[3 * i + 1];
              t[3 * i + 2] = vector.z * e[3 * i + 2];
            }

            return target;
          }
          /**
           * Solve Ax=b
           * @param b The right hand side
           * @param target Optional. Target vector to save in.
           * @return The solution x
           * @todo should reuse arrays
           */

          solve(b, target = new Vec3()) {
            // Construct equations
            const nr = 3; // num rows

            const nc = 4; // num cols

            const eqns = [];
            let i;
            let j;

            for (i = 0; i < nr * nc; i++) {
              eqns.push(0);
            }

            for (i = 0; i < 3; i++) {
              for (j = 0; j < 3; j++) {
                eqns[i + nc * j] = this.elements[i + 3 * j];
              }
            }

            eqns[3 + 4 * 0] = b.x;
            eqns[3 + 4 * 1] = b.y;
            eqns[3 + 4 * 2] = b.z; // Compute right upper triangular version of the matrix - Gauss elimination

            let n = 3;
            const k = n;
            let np;
            const kp = 4; // num rows

            let p;

            do {
              i = k - n;

              if (eqns[i + nc * i] === 0) {
                // the pivot is null, swap lines
                for (j = i + 1; j < k; j++) {
                  if (eqns[i + nc * j] !== 0) {
                    np = kp;

                    do {
                      // do ligne( i ) = ligne( i ) + ligne( k )
                      p = kp - np;
                      eqns[p + nc * i] += eqns[p + nc * j];
                    } while (--np);

                    break;
                  }
                }
              }

              if (eqns[i + nc * i] !== 0) {
                for (j = i + 1; j < k; j++) {
                  const multiplier = eqns[i + nc * j] / eqns[i + nc * i];
                  np = kp;

                  do {
                    // do ligne( k ) = ligne( k ) - multiplier * ligne( i )
                    p = kp - np;
                    eqns[p + nc * j] =
                      p <= i
                        ? 0
                        : eqns[p + nc * j] - eqns[p + nc * i] * multiplier;
                  } while (--np);
                }
              }
            } while (--n); // Get the solution

            target.z = eqns[2 * nc + 3] / eqns[2 * nc + 2];
            target.y =
              (eqns[1 * nc + 3] - eqns[1 * nc + 2] * target.z) /
              eqns[1 * nc + 1];
            target.x =
              (eqns[0 * nc + 3] -
                eqns[0 * nc + 2] * target.z -
                eqns[0 * nc + 1] * target.y) /
              eqns[0 * nc + 0];

            if (
              isNaN(target.x) ||
              isNaN(target.y) ||
              isNaN(target.z) ||
              target.x === Infinity ||
              target.y === Infinity ||
              target.z === Infinity
            ) {
              throw (
                "Could not solve equation! Got x=[" +
                target.toString() +
                "], b=[" +
                b.toString() +
                "], A=[" +
                this.toString() +
                "]"
              );
            }

            return target;
          }
          /**
           * Get an element in the matrix by index. Index starts at 0, not 1!!!
           * @param value If provided, the matrix element will be set to this value.
           */

          e(row, column, value) {
            if (value === undefined) {
              return this.elements[column + 3 * row];
            } else {
              // Set value
              this.elements[column + 3 * row] = value;
            }
          }
          /**
           * Copy another matrix into this matrix object.
           */

          copy(matrix) {
            for (let i = 0; i < matrix.elements.length; i++) {
              this.elements[i] = matrix.elements[i];
            }

            return this;
          }
          /**
           * Returns a string representation of the matrix.
           */

          toString() {
            let r = "";
            const sep = ",";

            for (let i = 0; i < 9; i++) {
              r += this.elements[i] + sep;
            }

            return r;
          }
          /**
           * reverse the matrix
           * @param target Target matrix to save in.
           * @return The solution x
           */

          reverse(target = new Mat3()) {
            // Construct equations
            const nr = 3; // num rows

            const nc = 6; // num cols

            const eqns = reverse_eqns;
            let i;
            let j;

            for (i = 0; i < 3; i++) {
              for (j = 0; j < 3; j++) {
                eqns[i + nc * j] = this.elements[i + 3 * j];
              }
            }

            eqns[3 + 6 * 0] = 1;
            eqns[3 + 6 * 1] = 0;
            eqns[3 + 6 * 2] = 0;
            eqns[4 + 6 * 0] = 0;
            eqns[4 + 6 * 1] = 1;
            eqns[4 + 6 * 2] = 0;
            eqns[5 + 6 * 0] = 0;
            eqns[5 + 6 * 1] = 0;
            eqns[5 + 6 * 2] = 1; // Compute right upper triangular version of the matrix - Gauss elimination

            let n = 3;
            const k = n;
            let np;
            const kp = nc; // num rows

            let p;

            do {
              i = k - n;

              if (eqns[i + nc * i] === 0) {
                // the pivot is null, swap lines
                for (j = i + 1; j < k; j++) {
                  if (eqns[i + nc * j] !== 0) {
                    np = kp;

                    do {
                      // do line( i ) = line( i ) + line( k )
                      p = kp - np;
                      eqns[p + nc * i] += eqns[p + nc * j];
                    } while (--np);

                    break;
                  }
                }
              }

              if (eqns[i + nc * i] !== 0) {
                for (j = i + 1; j < k; j++) {
                  const multiplier = eqns[i + nc * j] / eqns[i + nc * i];
                  np = kp;

                  do {
                    // do line( k ) = line( k ) - multiplier * line( i )
                    p = kp - np;
                    eqns[p + nc * j] =
                      p <= i
                        ? 0
                        : eqns[p + nc * j] - eqns[p + nc * i] * multiplier;
                  } while (--np);
                }
              }
            } while (--n); // eliminate the upper left triangle of the matrix

            i = 2;

            do {
              j = i - 1;

              do {
                const multiplier = eqns[i + nc * j] / eqns[i + nc * i];
                np = nc;

                do {
                  p = nc - np;
                  eqns[p + nc * j] =
                    eqns[p + nc * j] - eqns[p + nc * i] * multiplier;
                } while (--np);
              } while (j--);
            } while (--i); // operations on the diagonal

            i = 2;

            do {
              const multiplier = 1 / eqns[i + nc * i];
              np = nc;

              do {
                p = nc - np;
                eqns[p + nc * i] = eqns[p + nc * i] * multiplier;
              } while (--np);
            } while (i--);

            i = 2;

            do {
              j = 2;

              do {
                p = eqns[nr + j + nc * i];

                if (isNaN(p) || p === Infinity) {
                  throw "Could not reverse! A=[" + this.toString() + "]";
                }

                target.e(i, j, p);
              } while (j--);
            } while (i--);

            return target;
          }
          /**
           * Set the matrix from a quaterion
           */

          setRotationFromQuaternion(q) {
            const x = q.x;
            const y = q.y;
            const z = q.z;
            const w = q.w;
            const x2 = x + x;
            const y2 = y + y;
            const z2 = z + z;
            const xx = x * x2;
            const xy = x * y2;
            const xz = x * z2;
            const yy = y * y2;
            const yz = y * z2;
            const zz = z * z2;
            const wx = w * x2;
            const wy = w * y2;
            const wz = w * z2;
            const e = this.elements;
            e[3 * 0 + 0] = 1 - (yy + zz);
            e[3 * 0 + 1] = xy - wz;
            e[3 * 0 + 2] = xz + wy;
            e[3 * 1 + 0] = xy + wz;
            e[3 * 1 + 1] = 1 - (xx + zz);
            e[3 * 1 + 2] = yz - wx;
            e[3 * 2 + 0] = xz - wy;
            e[3 * 2 + 1] = yz + wx;
            e[3 * 2 + 2] = 1 - (xx + yy);
            return this;
          }
          /**
           * Transpose the matrix
           * @param target Optional. Where to store the result.
           * @return The target Mat3, or a new Mat3 if target was omitted.
           */

          transpose(target = new Mat3()) {
            const M = this.elements;
            const T = target.elements;
            let tmp; //Set diagonals

            T[0] = M[0];
            T[4] = M[4];
            T[8] = M[8];
            tmp = M[1];
            T[1] = M[3];
            T[3] = tmp;
            tmp = M[2];
            T[2] = M[6];
            T[6] = tmp;
            tmp = M[5];
            T[5] = M[7];
            T[7] = tmp;
            return target;
          }
        }
        const reverse_eqns = [
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];

        /**
         * 3-dimensional vector
         * @example
         *     const v = new Vec3(1, 2, 3)
         *     console.log('x=' + v.x) // x=1
         */

        class Vec3 {
          constructor(x = 0.0, y = 0.0, z = 0.0) {
            this.x = void 0;
            this.y = void 0;
            this.z = void 0;
            this.x = x;
            this.y = y;
            this.z = z;
          }
          /**
           * Vector cross product
           * @param target Optional target to save in.
           */

          cross(vector, target = new Vec3()) {
            const vx = vector.x;
            const vy = vector.y;
            const vz = vector.z;
            const x = this.x;
            const y = this.y;
            const z = this.z;
            target.x = y * vz - z * vy;
            target.y = z * vx - x * vz;
            target.z = x * vy - y * vx;
            return target;
          }
          /**
           * Set the vectors' 3 elements
           */

          set(x, y, z) {
            this.x = x;
            this.y = y;
            this.z = z;
            return this;
          }
          /**
           * Set all components of the vector to zero.
           */

          setZero() {
            this.x = this.y = this.z = 0;
          }
          /**
           * Vector addition
           */

          vadd(vector, target) {
            if (target) {
              target.x = vector.x + this.x;
              target.y = vector.y + this.y;
              target.z = vector.z + this.z;
            } else {
              return new Vec3(
                this.x + vector.x,
                this.y + vector.y,
                this.z + vector.z
              );
            }
          }
          /**
           * Vector subtraction
           * @param target Optional target to save in.
           */

          vsub(vector, target) {
            if (target) {
              target.x = this.x - vector.x;
              target.y = this.y - vector.y;
              target.z = this.z - vector.z;
            } else {
              return new Vec3(
                this.x - vector.x,
                this.y - vector.y,
                this.z - vector.z
              );
            }
          }
          /**
           * Get the cross product matrix a_cross from a vector, such that a x b = a_cross * b = c
           *
           * See {@link https://www8.cs.umu.se/kurser/TDBD24/VT06/lectures/Lecture6.pdf Umeå University Lecture}
           */

          crossmat() {
            return new Mat3([
              0,
              -this.z,
              this.y,
              this.z,
              0,
              -this.x,
              -this.y,
              this.x,
              0,
            ]);
          }
          /**
           * Normalize the vector. Note that this changes the values in the vector.
           * @return Returns the norm of the vector
           */

          normalize() {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            const n = Math.sqrt(x * x + y * y + z * z);

            if (n > 0.0) {
              const invN = 1 / n;
              this.x *= invN;
              this.y *= invN;
              this.z *= invN;
            } else {
              // Make something up
              this.x = 0;
              this.y = 0;
              this.z = 0;
            }

            return n;
          }
          /**
           * Get the version of this vector that is of length 1.
           * @param target Optional target to save in
           * @return Returns the unit vector
           */

          unit(target = new Vec3()) {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            let ninv = Math.sqrt(x * x + y * y + z * z);

            if (ninv > 0.0) {
              ninv = 1.0 / ninv;
              target.x = x * ninv;
              target.y = y * ninv;
              target.z = z * ninv;
            } else {
              target.x = 1;
              target.y = 0;
              target.z = 0;
            }

            return target;
          }
          /**
           * Get the length of the vector
           */

          length() {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            return Math.sqrt(x * x + y * y + z * z);
          }
          /**
           * Get the squared length of the vector.
           */

          lengthSquared() {
            return this.dot(this);
          }
          /**
           * Get distance from this point to another point
           */

          distanceTo(p) {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            const px = p.x;
            const py = p.y;
            const pz = p.z;
            return Math.sqrt(
              (px - x) * (px - x) + (py - y) * (py - y) + (pz - z) * (pz - z)
            );
          }
          /**
           * Get squared distance from this point to another point
           */

          distanceSquared(p) {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            const px = p.x;
            const py = p.y;
            const pz = p.z;
            return (
              (px - x) * (px - x) + (py - y) * (py - y) + (pz - z) * (pz - z)
            );
          }
          /**
           * Multiply all the components of the vector with a scalar.
           * @param target The vector to save the result in.
           */

          scale(scalar, target = new Vec3()) {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            target.x = scalar * x;
            target.y = scalar * y;
            target.z = scalar * z;
            return target;
          }
          /**
           * Multiply the vector with an other vector, component-wise.
           * @param target The vector to save the result in.
           */

          vmul(vector, target = new Vec3()) {
            target.x = vector.x * this.x;
            target.y = vector.y * this.y;
            target.z = vector.z * this.z;
            return target;
          }
          /**
           * Scale a vector and add it to this vector. Save the result in "target". (target = this + vector * scalar)
           * @param target The vector to save the result in.
           */

          addScaledVector(scalar, vector, target = new Vec3()) {
            target.x = this.x + scalar * vector.x;
            target.y = this.y + scalar * vector.y;
            target.z = this.z + scalar * vector.z;
            return target;
          }
          /**
           * Calculate dot product
           * @param vector
           */

          dot(vector) {
            return this.x * vector.x + this.y * vector.y + this.z * vector.z;
          }

          isZero() {
            return this.x === 0 && this.y === 0 && this.z === 0;
          }
          /**
           * Make the vector point in the opposite direction.
           * @param target Optional target to save in
           */

          negate(target = new Vec3()) {
            target.x = -this.x;
            target.y = -this.y;
            target.z = -this.z;
            return target;
          }
          /**
           * Compute two artificial tangents to the vector
           * @param t1 Vector object to save the first tangent in
           * @param t2 Vector object to save the second tangent in
           */

          tangents(t1, t2) {
            const norm = this.length();

            if (norm > 0.0) {
              const n = Vec3_tangents_n;
              const inorm = 1 / norm;
              n.set(this.x * inorm, this.y * inorm, this.z * inorm);
              const randVec = Vec3_tangents_randVec;

              if (Math.abs(n.x) < 0.9) {
                randVec.set(1, 0, 0);
                n.cross(randVec, t1);
              } else {
                randVec.set(0, 1, 0);
                n.cross(randVec, t1);
              }

              n.cross(t1, t2);
            } else {
              // The normal length is zero, make something up
              t1.set(1, 0, 0);
              t2.set(0, 1, 0);
            }
          }
          /**
           * Converts to a more readable format
           */

          toString() {
            return this.x + "," + this.y + "," + this.z;
          }
          /**
           * Converts to an array
           */

          toArray() {
            return [this.x, this.y, this.z];
          }
          /**
           * Copies value of source to this vector.
           */

          copy(vector) {
            this.x = vector.x;
            this.y = vector.y;
            this.z = vector.z;
            return this;
          }
          /**
           * Do a linear interpolation between two vectors
           * @param t A number between 0 and 1. 0 will make this function return u, and 1 will make it return v. Numbers in between will generate a vector in between them.
           */

          lerp(vector, t, target) {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            target.x = x + (vector.x - x) * t;
            target.y = y + (vector.y - y) * t;
            target.z = z + (vector.z - z) * t;
          }
          /**
           * Check if a vector equals is almost equal to another one.
           */

          almostEquals(vector, precision = 1e-6) {
            if (
              Math.abs(this.x - vector.x) > precision ||
              Math.abs(this.y - vector.y) > precision ||
              Math.abs(this.z - vector.z) > precision
            ) {
              return false;
            }

            return true;
          }
          /**
           * Check if a vector is almost zero
           */

          almostZero(precision = 1e-6) {
            if (
              Math.abs(this.x) > precision ||
              Math.abs(this.y) > precision ||
              Math.abs(this.z) > precision
            ) {
              return false;
            }

            return true;
          }
          /**
           * Check if the vector is anti-parallel to another vector.
           * @param precision Set to zero for exact comparisons
           */

          isAntiparallelTo(vector, precision) {
            this.negate(antip_neg);
            return antip_neg.almostEquals(vector, precision);
          }
          /**
           * Clone the vector
           */

          clone() {
            return new Vec3(this.x, this.y, this.z);
          }
        }
        Vec3.ZERO = void 0;
        Vec3.UNIT_X = void 0;
        Vec3.UNIT_Y = void 0;
        Vec3.UNIT_Z = void 0;
        Vec3.ZERO = new Vec3(0, 0, 0);
        Vec3.UNIT_X = new Vec3(1, 0, 0);
        Vec3.UNIT_Y = new Vec3(0, 1, 0);
        Vec3.UNIT_Z = new Vec3(0, 0, 1);
        const Vec3_tangents_n = new Vec3();
        const Vec3_tangents_randVec = new Vec3();
        const antip_neg = new Vec3();

        /**
         * Axis aligned bounding box class.
         */
        class AABB {
          /**
           * The lower bound of the bounding box
           */

          /**
           * The upper bound of the bounding box
           */
          constructor(options = {}) {
            this.lowerBound = void 0;
            this.upperBound = void 0;
            this.lowerBound = new Vec3();
            this.upperBound = new Vec3();

            if (options.lowerBound) {
              this.lowerBound.copy(options.lowerBound);
            }

            if (options.upperBound) {
              this.upperBound.copy(options.upperBound);
            }
          }
          /**
           * Set the AABB bounds from a set of points.
           * @param points An array of Vec3's.
           * @return The self object
           */

          setFromPoints(points, position, quaternion, skinSize) {
            const l = this.lowerBound;
            const u = this.upperBound;
            const q = quaternion; // Set to the first point

            l.copy(points[0]);

            if (q) {
              q.vmult(l, l);
            }

            u.copy(l);

            for (let i = 1; i < points.length; i++) {
              let p = points[i];

              if (q) {
                q.vmult(p, tmp$1);
                p = tmp$1;
              }

              if (p.x > u.x) {
                u.x = p.x;
              }

              if (p.x < l.x) {
                l.x = p.x;
              }

              if (p.y > u.y) {
                u.y = p.y;
              }

              if (p.y < l.y) {
                l.y = p.y;
              }

              if (p.z > u.z) {
                u.z = p.z;
              }

              if (p.z < l.z) {
                l.z = p.z;
              }
            } // Add offset

            if (position) {
              position.vadd(l, l);
              position.vadd(u, u);
            }

            if (skinSize) {
              l.x -= skinSize;
              l.y -= skinSize;
              l.z -= skinSize;
              u.x += skinSize;
              u.y += skinSize;
              u.z += skinSize;
            }

            return this;
          }
          /**
           * Copy bounds from an AABB to this AABB
           * @param aabb Source to copy from
           * @return The this object, for chainability
           */

          copy(aabb) {
            this.lowerBound.copy(aabb.lowerBound);
            this.upperBound.copy(aabb.upperBound);
            return this;
          }
          /**
           * Clone an AABB
           */

          clone() {
            return new AABB().copy(this);
          }
          /**
           * Extend this AABB so that it covers the given AABB too.
           */

          extend(aabb) {
            this.lowerBound.x = Math.min(this.lowerBound.x, aabb.lowerBound.x);
            this.upperBound.x = Math.max(this.upperBound.x, aabb.upperBound.x);
            this.lowerBound.y = Math.min(this.lowerBound.y, aabb.lowerBound.y);
            this.upperBound.y = Math.max(this.upperBound.y, aabb.upperBound.y);
            this.lowerBound.z = Math.min(this.lowerBound.z, aabb.lowerBound.z);
            this.upperBound.z = Math.max(this.upperBound.z, aabb.upperBound.z);
          }
          /**
           * Returns true if the given AABB overlaps this AABB.
           */

          overlaps(aabb) {
            const l1 = this.lowerBound;
            const u1 = this.upperBound;
            const l2 = aabb.lowerBound;
            const u2 = aabb.upperBound; //      l2        u2
            //      |---------|
            // |--------|
            // l1       u1

            const overlapsX =
              (l2.x <= u1.x && u1.x <= u2.x) || (l1.x <= u2.x && u2.x <= u1.x);
            const overlapsY =
              (l2.y <= u1.y && u1.y <= u2.y) || (l1.y <= u2.y && u2.y <= u1.y);
            const overlapsZ =
              (l2.z <= u1.z && u1.z <= u2.z) || (l1.z <= u2.z && u2.z <= u1.z);
            return overlapsX && overlapsY && overlapsZ;
          } // Mostly for debugging

          volume() {
            const l = this.lowerBound;
            const u = this.upperBound;
            return (u.x - l.x) * (u.y - l.y) * (u.z - l.z);
          }
          /**
           * Returns true if the given AABB is fully contained in this AABB.
           */

          contains(aabb) {
            const l1 = this.lowerBound;
            const u1 = this.upperBound;
            const l2 = aabb.lowerBound;
            const u2 = aabb.upperBound; //      l2        u2
            //      |---------|
            // |---------------|
            // l1              u1

            return (
              l1.x <= l2.x &&
              u1.x >= u2.x &&
              l1.y <= l2.y &&
              u1.y >= u2.y &&
              l1.z <= l2.z &&
              u1.z >= u2.z
            );
          }

          getCorners(a, b, c, d, e, f, g, h) {
            const l = this.lowerBound;
            const u = this.upperBound;
            a.copy(l);
            b.set(u.x, l.y, l.z);
            c.set(u.x, u.y, l.z);
            d.set(l.x, u.y, u.z);
            e.set(u.x, l.y, u.z);
            f.set(l.x, u.y, l.z);
            g.set(l.x, l.y, u.z);
            h.copy(u);
          }
          /**
           * Get the representation of an AABB in another frame.
           * @return The "target" AABB object.
           */

          toLocalFrame(frame, target) {
            const corners = transformIntoFrame_corners;
            const a = corners[0];
            const b = corners[1];
            const c = corners[2];
            const d = corners[3];
            const e = corners[4];
            const f = corners[5];
            const g = corners[6];
            const h = corners[7]; // Get corners in current frame

            this.getCorners(a, b, c, d, e, f, g, h); // Transform them to new local frame

            for (let i = 0; i !== 8; i++) {
              const corner = corners[i];
              frame.pointToLocal(corner, corner);
            }

            return target.setFromPoints(corners);
          }
          /**
           * Get the representation of an AABB in the global frame.
           * @return The "target" AABB object.
           */

          toWorldFrame(frame, target) {
            const corners = transformIntoFrame_corners;
            const a = corners[0];
            const b = corners[1];
            const c = corners[2];
            const d = corners[3];
            const e = corners[4];
            const f = corners[5];
            const g = corners[6];
            const h = corners[7]; // Get corners in current frame

            this.getCorners(a, b, c, d, e, f, g, h); // Transform them to new local frame

            for (let i = 0; i !== 8; i++) {
              const corner = corners[i];
              frame.pointToWorld(corner, corner);
            }

            return target.setFromPoints(corners);
          }
          /**
           * Check if the AABB is hit by a ray.
           */

          overlapsRay(ray) {
            const { direction, from } = ray; // const t = 0
            // ray.direction is unit direction vector of ray

            const dirFracX = 1 / direction.x;
            const dirFracY = 1 / direction.y;
            const dirFracZ = 1 / direction.z; // this.lowerBound is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner

            const t1 = (this.lowerBound.x - from.x) * dirFracX;
            const t2 = (this.upperBound.x - from.x) * dirFracX;
            const t3 = (this.lowerBound.y - from.y) * dirFracY;
            const t4 = (this.upperBound.y - from.y) * dirFracY;
            const t5 = (this.lowerBound.z - from.z) * dirFracZ;
            const t6 = (this.upperBound.z - from.z) * dirFracZ; // const tmin = Math.max(Math.max(Math.min(t1, t2), Math.min(t3, t4)));
            // const tmax = Math.min(Math.min(Math.max(t1, t2), Math.max(t3, t4)));

            const tmin = Math.max(
              Math.max(Math.min(t1, t2), Math.min(t3, t4)),
              Math.min(t5, t6)
            );
            const tmax = Math.min(
              Math.min(Math.max(t1, t2), Math.max(t3, t4)),
              Math.max(t5, t6)
            ); // if tmax < 0, ray (line) is intersecting AABB, but whole AABB is behing us

            if (tmax < 0) {
              //t = tmax;
              return false;
            } // if tmin > tmax, ray doesn't intersect AABB

            if (tmin > tmax) {
              //t = tmax;
              return false;
            }

            return true;
          }
        }
        const tmp$1 = new Vec3();
        const transformIntoFrame_corners = [
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
        ];

        /**
         * Collision "matrix".
         * It's actually a triangular-shaped array of whether two bodies are touching this step, for reference next step
         */
        class ArrayCollisionMatrix {
          /**
           * The matrix storage.
           */
          constructor() {
            this.matrix = void 0;
            this.matrix = [];
          }
          /**
           * Get an element
           */

          get(bi, bj) {
            let { index: i } = bi;
            let { index: j } = bj;

            if (j > i) {
              const temp = j;
              j = i;
              i = temp;
            }

            return this.matrix[((i * (i + 1)) >> 1) + j - 1];
          }
          /**
           * Set an element
           */

          set(bi, bj, value) {
            let { index: i } = bi;
            let { index: j } = bj;

            if (j > i) {
              const temp = j;
              j = i;
              i = temp;
            }

            this.matrix[((i * (i + 1)) >> 1) + j - 1] = value ? 1 : 0;
          }
          /**
           * Sets all elements to zero
           */

          reset() {
            for (let i = 0, l = this.matrix.length; i !== l; i++) {
              this.matrix[i] = 0;
            }
          }
          /**
           * Sets the max number of objects
           */

          setNumObjects(n) {
            this.matrix.length = (n * (n - 1)) >> 1;
          }
        }

        /**
         * Base class for objects that dispatches events.
         */
        class EventTarget {
          constructor() {
            this._listeners = void 0;
          }

          /**
           * Add an event listener
           * @return The self object, for chainability.
           */
          addEventListener(type, listener) {
            if (this._listeners === undefined) {
              this._listeners = {};
            }

            const listeners = this._listeners;

            if (listeners[type] === undefined) {
              listeners[type] = [];
            }

            if (!listeners[type].includes(listener)) {
              listeners[type].push(listener);
            }

            return this;
          }
          /**
           * Check if an event listener is added
           */

          hasEventListener(type, listener) {
            if (this._listeners === undefined) {
              return false;
            }

            const listeners = this._listeners;

            if (
              listeners[type] !== undefined &&
              listeners[type].includes(listener)
            ) {
              return true;
            }

            return false;
          }
          /**
           * Check if any event listener of the given type is added
           */

          hasAnyEventListener(type) {
            if (this._listeners === undefined) {
              return false;
            }

            const listeners = this._listeners;
            return listeners[type] !== undefined;
          }
          /**
           * Remove an event listener
           * @return The self object, for chainability.
           */

          removeEventListener(type, listener) {
            if (this._listeners === undefined) {
              return this;
            }

            const listeners = this._listeners;

            if (listeners[type] === undefined) {
              return this;
            }

            const index = listeners[type].indexOf(listener);

            if (index !== -1) {
              listeners[type].splice(index, 1);
            }

            return this;
          }
          /**
           * Emit an event.
           * @return The self object, for chainability.
           */

          dispatchEvent(event) {
            if (this._listeners === undefined) {
              return this;
            }

            const listeners = this._listeners;
            const listenerArray = listeners[event.type];

            if (listenerArray !== undefined) {
              event.target = this;

              for (let i = 0, l = listenerArray.length; i < l; i++) {
                listenerArray[i].call(this, event);
              }
            }

            return this;
          }
        }

        /**
         * A Quaternion describes a rotation in 3D space. The Quaternion is mathematically defined as Q = x*i + y*j + z*k + w, where (i,j,k) are imaginary basis vectors. (x,y,z) can be seen as a vector related to the axis of rotation, while the real multiplier, w, is related to the amount of rotation.
         * @param x Multiplier of the imaginary basis vector i.
         * @param y Multiplier of the imaginary basis vector j.
         * @param z Multiplier of the imaginary basis vector k.
         * @param w Multiplier of the real part.
         * @see http://en.wikipedia.org/wiki/Quaternion
         */

        class Quaternion {
          constructor(x = 0, y = 0, z = 0, w = 1) {
            this.x = void 0;
            this.y = void 0;
            this.z = void 0;
            this.w = void 0;
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
          }
          /**
           * Set the value of the quaternion.
           */

          set(x, y, z, w) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
            return this;
          }
          /**
           * Convert to a readable format
           * @return "x,y,z,w"
           */

          toString() {
            return this.x + "," + this.y + "," + this.z + "," + this.w;
          }
          /**
           * Convert to an Array
           * @return [x, y, z, w]
           */

          toArray() {
            return [this.x, this.y, this.z, this.w];
          }
          /**
           * Set the quaternion components given an axis and an angle in radians.
           */

          setFromAxisAngle(vector, angle) {
            const s = Math.sin(angle * 0.5);
            this.x = vector.x * s;
            this.y = vector.y * s;
            this.z = vector.z * s;
            this.w = Math.cos(angle * 0.5);
            return this;
          }
          /**
           * Converts the quaternion to [ axis, angle ] representation.
           * @param targetAxis A vector object to reuse for storing the axis.
           * @return An array, first element is the axis and the second is the angle in radians.
           */

          toAxisAngle(targetAxis = new Vec3()) {
            this.normalize(); // if w>1 acos and sqrt will produce errors, this cant happen if quaternion is normalised

            const angle = 2 * Math.acos(this.w);
            const s = Math.sqrt(1 - this.w * this.w); // assuming quaternion normalised then w is less than 1, so term always positive.

            if (s < 0.001) {
              // test to avoid divide by zero, s is always positive due to sqrt
              // if s close to zero then direction of axis not important
              targetAxis.x = this.x; // if it is important that axis is normalised then replace with x=1; y=z=0;

              targetAxis.y = this.y;
              targetAxis.z = this.z;
            } else {
              targetAxis.x = this.x / s; // normalise axis

              targetAxis.y = this.y / s;
              targetAxis.z = this.z / s;
            }

            return [targetAxis, angle];
          }
          /**
           * Set the quaternion value given two vectors. The resulting rotation will be the needed rotation to rotate u to v.
           */

          setFromVectors(u, v) {
            if (u.isAntiparallelTo(v)) {
              const t1 = sfv_t1;
              const t2 = sfv_t2;
              u.tangents(t1, t2);
              this.setFromAxisAngle(t1, Math.PI);
            } else {
              const a = u.cross(v);
              this.x = a.x;
              this.y = a.y;
              this.z = a.z;
              this.w = Math.sqrt(u.length() ** 2 * v.length() ** 2) + u.dot(v);
              this.normalize();
            }

            return this;
          }
          /**
           * Multiply the quaternion with an other quaternion.
           */

          mult(quat, target = new Quaternion()) {
            const ax = this.x;
            const ay = this.y;
            const az = this.z;
            const aw = this.w;
            const bx = quat.x;
            const by = quat.y;
            const bz = quat.z;
            const bw = quat.w;
            target.x = ax * bw + aw * bx + ay * bz - az * by;
            target.y = ay * bw + aw * by + az * bx - ax * bz;
            target.z = az * bw + aw * bz + ax * by - ay * bx;
            target.w = aw * bw - ax * bx - ay * by - az * bz;
            return target;
          }
          /**
           * Get the inverse quaternion rotation.
           */

          inverse(target = new Quaternion()) {
            const x = this.x;
            const y = this.y;
            const z = this.z;
            const w = this.w;
            this.conjugate(target);
            const inorm2 = 1 / (x * x + y * y + z * z + w * w);
            target.x *= inorm2;
            target.y *= inorm2;
            target.z *= inorm2;
            target.w *= inorm2;
            return target;
          }
          /**
           * Get the quaternion conjugate
           */

          conjugate(target = new Quaternion()) {
            target.x = -this.x;
            target.y = -this.y;
            target.z = -this.z;
            target.w = this.w;
            return target;
          }
          /**
           * Normalize the quaternion. Note that this changes the values of the quaternion.
           */

          normalize() {
            let l = Math.sqrt(
              this.x * this.x +
                this.y * this.y +
                this.z * this.z +
                this.w * this.w
            );

            if (l === 0) {
              this.x = 0;
              this.y = 0;
              this.z = 0;
              this.w = 0;
            } else {
              l = 1 / l;
              this.x *= l;
              this.y *= l;
              this.z *= l;
              this.w *= l;
            }

            return this;
          }
          /**
           * Approximation of quaternion normalization. Works best when quat is already almost-normalized.
           * @author unphased, https://github.com/unphased
           */

          normalizeFast() {
            const f =
              (3.0 -
                (this.x * this.x +
                  this.y * this.y +
                  this.z * this.z +
                  this.w * this.w)) /
              2.0;

            if (f === 0) {
              this.x = 0;
              this.y = 0;
              this.z = 0;
              this.w = 0;
            } else {
              this.x *= f;
              this.y *= f;
              this.z *= f;
              this.w *= f;
            }

            return this;
          }
          /**
           * Multiply the quaternion by a vector
           */

          vmult(v, target = new Vec3()) {
            const x = v.x;
            const y = v.y;
            const z = v.z;
            const qx = this.x;
            const qy = this.y;
            const qz = this.z;
            const qw = this.w; // q*v

            const ix = qw * x + qy * z - qz * y;
            const iy = qw * y + qz * x - qx * z;
            const iz = qw * z + qx * y - qy * x;
            const iw = -qx * x - qy * y - qz * z;
            target.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
            target.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
            target.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;
            return target;
          }
          /**
           * Copies value of source to this quaternion.
           * @return this
           */

          copy(quat) {
            this.x = quat.x;
            this.y = quat.y;
            this.z = quat.z;
            this.w = quat.w;
            return this;
          }
          /**
           * Convert the quaternion to euler angle representation. Order: YZX, as this page describes: https://www.euclideanspace.com/maths/standards/index.htm
           * @param order Three-character string, defaults to "YZX"
           */

          toEuler(target, order = "YZX") {
            let heading;
            let attitude;
            let bank;
            const x = this.x;
            const y = this.y;
            const z = this.z;
            const w = this.w;

            switch (order) {
              case "YZX":
                const test = x * y + z * w;

                if (test > 0.499) {
                  // singularity at north pole
                  heading = 2 * Math.atan2(x, w);
                  attitude = Math.PI / 2;
                  bank = 0;
                }

                if (test < -0.499) {
                  // singularity at south pole
                  heading = -2 * Math.atan2(x, w);
                  attitude = -Math.PI / 2;
                  bank = 0;
                }

                if (heading === undefined) {
                  const sqx = x * x;
                  const sqy = y * y;
                  const sqz = z * z;
                  heading = Math.atan2(
                    2 * y * w - 2 * x * z,
                    1 - 2 * sqy - 2 * sqz
                  ); // Heading

                  attitude = Math.asin(2 * test); // attitude

                  bank = Math.atan2(
                    2 * x * w - 2 * y * z,
                    1 - 2 * sqx - 2 * sqz
                  ); // bank
                }

                break;

              default:
                throw new Error("Euler order " + order + " not supported yet.");
            }

            target.y = heading;
            target.z = attitude;
            target.x = bank;
          }
          /**
           * @param order The order to apply angles: 'XYZ' or 'YXZ' or any other combination.
           *
           * See {@link https://www.mathworks.com/matlabcentral/fileexchange/20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors MathWorks} reference
           */

          setFromEuler(x, y, z, order = "XYZ") {
            const c1 = Math.cos(x / 2);
            const c2 = Math.cos(y / 2);
            const c3 = Math.cos(z / 2);
            const s1 = Math.sin(x / 2);
            const s2 = Math.sin(y / 2);
            const s3 = Math.sin(z / 2);

            if (order === "XYZ") {
              this.x = s1 * c2 * c3 + c1 * s2 * s3;
              this.y = c1 * s2 * c3 - s1 * c2 * s3;
              this.z = c1 * c2 * s3 + s1 * s2 * c3;
              this.w = c1 * c2 * c3 - s1 * s2 * s3;
            } else if (order === "YXZ") {
              this.x = s1 * c2 * c3 + c1 * s2 * s3;
              this.y = c1 * s2 * c3 - s1 * c2 * s3;
              this.z = c1 * c2 * s3 - s1 * s2 * c3;
              this.w = c1 * c2 * c3 + s1 * s2 * s3;
            } else if (order === "ZXY") {
              this.x = s1 * c2 * c3 - c1 * s2 * s3;
              this.y = c1 * s2 * c3 + s1 * c2 * s3;
              this.z = c1 * c2 * s3 + s1 * s2 * c3;
              this.w = c1 * c2 * c3 - s1 * s2 * s3;
            } else if (order === "ZYX") {
              this.x = s1 * c2 * c3 - c1 * s2 * s3;
              this.y = c1 * s2 * c3 + s1 * c2 * s3;
              this.z = c1 * c2 * s3 - s1 * s2 * c3;
              this.w = c1 * c2 * c3 + s1 * s2 * s3;
            } else if (order === "YZX") {
              this.x = s1 * c2 * c3 + c1 * s2 * s3;
              this.y = c1 * s2 * c3 + s1 * c2 * s3;
              this.z = c1 * c2 * s3 - s1 * s2 * c3;
              this.w = c1 * c2 * c3 - s1 * s2 * s3;
            } else if (order === "XZY") {
              this.x = s1 * c2 * c3 - c1 * s2 * s3;
              this.y = c1 * s2 * c3 - s1 * c2 * s3;
              this.z = c1 * c2 * s3 + s1 * s2 * c3;
              this.w = c1 * c2 * c3 + s1 * s2 * s3;
            }

            return this;
          }

          clone() {
            return new Quaternion(this.x, this.y, this.z, this.w);
          }
          /**
           * Performs a spherical linear interpolation between two quat
           *
           * @param toQuat second operand
           * @param t interpolation amount between the self quaternion and toQuat
           * @param target A quaternion to store the result in. If not provided, a new one will be created.
           * @returns {Quaternion} The "target" object
           */

          slerp(toQuat, t, target = new Quaternion()) {
            const ax = this.x;
            const ay = this.y;
            const az = this.z;
            const aw = this.w;
            let bx = toQuat.x;
            let by = toQuat.y;
            let bz = toQuat.z;
            let bw = toQuat.w;
            let omega;
            let cosom;
            let sinom;
            let scale0;
            let scale1; // calc cosine

            cosom = ax * bx + ay * by + az * bz + aw * bw; // adjust signs (if necessary)

            if (cosom < 0.0) {
              cosom = -cosom;
              bx = -bx;
              by = -by;
              bz = -bz;
              bw = -bw;
            } // calculate coefficients

            if (1.0 - cosom > 0.000001) {
              // standard case (slerp)
              omega = Math.acos(cosom);
              sinom = Math.sin(omega);
              scale0 = Math.sin((1.0 - t) * omega) / sinom;
              scale1 = Math.sin(t * omega) / sinom;
            } else {
              // "from" and "to" quaternions are very close
              //  ... so we can do a linear interpolation
              scale0 = 1.0 - t;
              scale1 = t;
            } // calculate final values

            target.x = scale0 * ax + scale1 * bx;
            target.y = scale0 * ay + scale1 * by;
            target.z = scale0 * az + scale1 * bz;
            target.w = scale0 * aw + scale1 * bw;
            return target;
          }
          /**
           * Rotate an absolute orientation quaternion given an angular velocity and a time step.
           */

          integrate(
            angularVelocity,
            dt,
            angularFactor,
            target = new Quaternion()
          ) {
            const ax = angularVelocity.x * angularFactor.x,
              ay = angularVelocity.y * angularFactor.y,
              az = angularVelocity.z * angularFactor.z,
              bx = this.x,
              by = this.y,
              bz = this.z,
              bw = this.w;
            const half_dt = dt * 0.5;
            target.x += half_dt * (ax * bw + ay * bz - az * by);
            target.y += half_dt * (ay * bw + az * bx - ax * bz);
            target.z += half_dt * (az * bw + ax * by - ay * bx);
            target.w += half_dt * (-ax * bx - ay * by - az * bz);
            return target;
          }
        }
        const sfv_t1 = new Vec3();
        const sfv_t2 = new Vec3();

        /**
         * The available shape types.
         */
        const SHAPE_TYPES = {
          /** SPHERE */
          SPHERE: 1,

          /** PLANE */
          PLANE: 2,

          /** BOX */
          BOX: 4,

          /** COMPOUND */
          COMPOUND: 8,

          /** CONVEXPOLYHEDRON */
          CONVEXPOLYHEDRON: 16,

          /** HEIGHTFIELD */
          HEIGHTFIELD: 32,

          /** PARTICLE */
          PARTICLE: 64,

          /** CYLINDER */
          CYLINDER: 128,

          /** TRIMESH */
          TRIMESH: 256,
        };
        /**
         * ShapeType
         */

        /**
         * Base class for shapes
         */
        class Shape {
          /**
           * Identifier of the Shape.
           */

          /**
           * The type of this shape. Must be set to an int > 0 by subclasses.
           */

          /**
           * The local bounding sphere radius of this shape.
           */

          /**
           * Whether to produce contact forces when in contact with other bodies. Note that contacts will be generated, but they will be disabled.
           * @default true
           */

          /**
           * @default 1
           */

          /**
           * @default -1
           */

          /**
           * Optional material of the shape that regulates contact properties.
           */

          /**
           * The body to which the shape is added to.
           */

          /**
           * All the Shape types.
           */
          constructor(options = {}) {
            this.id = void 0;
            this.type = void 0;
            this.boundingSphereRadius = void 0;
            this.collisionResponse = void 0;
            this.collisionFilterGroup = void 0;
            this.collisionFilterMask = void 0;
            this.material = void 0;
            this.body = void 0;
            this.id = Shape.idCounter++;
            this.type = options.type || 0;
            this.boundingSphereRadius = 0;
            this.collisionResponse = options.collisionResponse
              ? options.collisionResponse
              : true;
            this.collisionFilterGroup =
              options.collisionFilterGroup !== undefined
                ? options.collisionFilterGroup
                : 1;
            this.collisionFilterMask =
              options.collisionFilterMask !== undefined
                ? options.collisionFilterMask
                : -1;
            this.material = options.material ? options.material : null;
            this.body = null;
          }
          /**
           * Computes the bounding sphere radius.
           * The result is stored in the property `.boundingSphereRadius`
           */

          updateBoundingSphereRadius() {
            throw (
              "computeBoundingSphereRadius() not implemented for shape type " +
              this.type
            );
          }
          /**
           * Get the volume of this shape
           */

          volume() {
            throw "volume() not implemented for shape type " + this.type;
          }
          /**
           * Calculates the inertia in the local frame for this shape.
           * @see http://en.wikipedia.org/wiki/List_of_moments_of_inertia
           */

          calculateLocalInertia(mass, target) {
            throw (
              "calculateLocalInertia() not implemented for shape type " +
              this.type
            );
          }
          /**
           * @todo use abstract for these kind of methods
           */

          calculateWorldAABB(pos, quat, min, max) {
            throw (
              "calculateWorldAABB() not implemented for shape type " + this.type
            );
          }
        }
        Shape.idCounter = 0;
        Shape.types = SHAPE_TYPES;

        /**
         * Transformation utilities.
         */
        class Transform {
          /**
           * position
           */

          /**
           * quaternion
           */
          constructor(options = {}) {
            this.position = void 0;
            this.quaternion = void 0;
            this.position = new Vec3();
            this.quaternion = new Quaternion();

            if (options.position) {
              this.position.copy(options.position);
            }

            if (options.quaternion) {
              this.quaternion.copy(options.quaternion);
            }
          }
          /**
           * Get a global point in local transform coordinates.
           */

          pointToLocal(worldPoint, result) {
            return Transform.pointToLocalFrame(
              this.position,
              this.quaternion,
              worldPoint,
              result
            );
          }
          /**
           * Get a local point in global transform coordinates.
           */

          pointToWorld(localPoint, result) {
            return Transform.pointToWorldFrame(
              this.position,
              this.quaternion,
              localPoint,
              result
            );
          }
          /**
           * vectorToWorldFrame
           */

          vectorToWorldFrame(localVector, result = new Vec3()) {
            this.quaternion.vmult(localVector, result);
            return result;
          }
          /**
           * pointToLocalFrame
           */

          static pointToLocalFrame(
            position,
            quaternion,
            worldPoint,
            result = new Vec3()
          ) {
            worldPoint.vsub(position, result);
            quaternion.conjugate(tmpQuat$1);
            tmpQuat$1.vmult(result, result);
            return result;
          }
          /**
           * pointToWorldFrame
           */

          static pointToWorldFrame(
            position,
            quaternion,
            localPoint,
            result = new Vec3()
          ) {
            quaternion.vmult(localPoint, result);
            result.vadd(position, result);
            return result;
          }
          /**
           * vectorToWorldFrame
           */

          static vectorToWorldFrame(
            quaternion,
            localVector,
            result = new Vec3()
          ) {
            quaternion.vmult(localVector, result);
            return result;
          }
          /**
           * vectorToLocalFrame
           */

          static vectorToLocalFrame(
            position,
            quaternion,
            worldVector,
            result = new Vec3()
          ) {
            quaternion.w *= -1;
            quaternion.vmult(worldVector, result);
            quaternion.w *= -1;
            return result;
          }
        }
        const tmpQuat$1 = new Quaternion();

        /**
         * A set of polygons describing a convex shape.
         *
         * The shape MUST be convex for the code to work properly. No polygons may be coplanar (contained
         * in the same 3D plane), instead these should be merged into one polygon.
         *
         * @author qiao / https://github.com/qiao (original author, see https://github.com/qiao/three.js/commit/85026f0c769e4000148a67d45a9e9b9c5108836f)
         * @author schteppe / https://github.com/schteppe
         * @see https://www.altdevblogaday.com/2011/05/13/contact-generation-between-3d-convex-meshes/
         *
         * @todo Move the clipping functions to ContactGenerator?
         * @todo Automatically merge coplanar polygons in constructor.
         * @example
         *     const convexShape = new CANNON.ConvexPolyhedron({ vertices, faces })
         *     const convexBody = new CANNON.Body({ mass: 1, shape: convexShape })
         *     world.addBody(convexBody)
         */
        class ConvexPolyhedron extends Shape {
          /** vertices */

          /**
           * Array of integer arrays, indicating which vertices each face consists of
           */

          /** faceNormals */

          /** worldVertices */

          /** worldVerticesNeedsUpdate */

          /** worldFaceNormals */

          /** worldFaceNormalsNeedsUpdate */

          /**
           * If given, these locally defined, normalized axes are the only ones being checked when doing separating axis check.
           */

          /** uniqueEdges */

          /**
           * @param vertices An array of Vec3's
           * @param faces Array of integer arrays, describing which vertices that is included in each face.
           */
          constructor(props = {}) {
            const {
              vertices = [],
              faces = [],
              normals = [],
              axes,
              boundingSphereRadius,
            } = props;
            super({
              type: Shape.types.CONVEXPOLYHEDRON,
            });
            this.vertices = void 0;
            this.faces = void 0;
            this.faceNormals = void 0;
            this.worldVertices = void 0;
            this.worldVerticesNeedsUpdate = void 0;
            this.worldFaceNormals = void 0;
            this.worldFaceNormalsNeedsUpdate = void 0;
            this.uniqueAxes = void 0;
            this.uniqueEdges = void 0;
            this.vertices = vertices;
            this.faces = faces;
            this.faceNormals = normals;

            if (this.faceNormals.length === 0) {
              this.computeNormals();
            }

            if (!boundingSphereRadius) {
              this.updateBoundingSphereRadius();
            } else {
              this.boundingSphereRadius = boundingSphereRadius;
            }

            this.worldVertices = []; // World transformed version of .vertices

            this.worldVerticesNeedsUpdate = true;
            this.worldFaceNormals = []; // World transformed version of .faceNormals

            this.worldFaceNormalsNeedsUpdate = true;
            this.uniqueAxes = axes ? axes.slice() : null;
            this.uniqueEdges = [];
            this.computeEdges();
          }
          /**
           * Computes uniqueEdges
           */

          computeEdges() {
            const faces = this.faces;
            const vertices = this.vertices;
            const edges = this.uniqueEdges;
            edges.length = 0;
            const edge = new Vec3();

            for (let i = 0; i !== faces.length; i++) {
              const face = faces[i];
              const numVertices = face.length;

              for (let j = 0; j !== numVertices; j++) {
                const k = (j + 1) % numVertices;
                vertices[face[j]].vsub(vertices[face[k]], edge);
                edge.normalize();
                let found = false;

                for (let p = 0; p !== edges.length; p++) {
                  if (
                    edges[p].almostEquals(edge) ||
                    edges[p].almostEquals(edge)
                  ) {
                    found = true;
                    break;
                  }
                }

                if (!found) {
                  edges.push(edge.clone());
                }
              }
            }
          }
          /**
           * Compute the normals of the faces.
           * Will reuse existing Vec3 objects in the `faceNormals` array if they exist.
           */

          computeNormals() {
            this.faceNormals.length = this.faces.length; // Generate normals

            for (let i = 0; i < this.faces.length; i++) {
              // Check so all vertices exists for this face
              for (let j = 0; j < this.faces[i].length; j++) {
                if (!this.vertices[this.faces[i][j]]) {
                  throw new Error("Vertex " + this.faces[i][j] + " not found!");
                }
              }

              const n = this.faceNormals[i] || new Vec3();
              this.getFaceNormal(i, n);
              n.negate(n);
              this.faceNormals[i] = n;
              const vertex = this.vertices[this.faces[i][0]];

              if (n.dot(vertex) < 0) {
                console.error(
                  ".faceNormals[" +
                    i +
                    "] = Vec3(" +
                    n.toString() +
                    ") looks like it points into the shape? The vertices follow. Make sure they are ordered CCW around the normal, using the right hand rule."
                );

                for (let j = 0; j < this.faces[i].length; j++) {
                  console.warn(
                    ".vertices[" +
                      this.faces[i][j] +
                      "] = Vec3(" +
                      this.vertices[this.faces[i][j]].toString() +
                      ")"
                  );
                }
              }
            }
          }
          /**
           * Compute the normal of a face from its vertices
           */

          getFaceNormal(i, target) {
            const f = this.faces[i];
            const va = this.vertices[f[0]];
            const vb = this.vertices[f[1]];
            const vc = this.vertices[f[2]];
            ConvexPolyhedron.computeNormal(va, vb, vc, target);
          }
          /**
           * Get face normal given 3 vertices
           */

          static computeNormal(va, vb, vc, target) {
            const cb = new Vec3();
            const ab = new Vec3();
            vb.vsub(va, ab);
            vc.vsub(vb, cb);
            cb.cross(ab, target);

            if (!target.isZero()) {
              target.normalize();
            }
          }
          /**
           * @param minDist Clamp distance
           * @param result The an array of contact point objects, see clipFaceAgainstHull
           */

          clipAgainstHull(
            posA,
            quatA,
            hullB,
            posB,
            quatB,
            separatingNormal,
            minDist,
            maxDist,
            result
          ) {
            const WorldNormal = new Vec3();
            let closestFaceB = -1;
            let dmax = -Number.MAX_VALUE;

            for (let face = 0; face < hullB.faces.length; face++) {
              WorldNormal.copy(hullB.faceNormals[face]);
              quatB.vmult(WorldNormal, WorldNormal);
              const d = WorldNormal.dot(separatingNormal);

              if (d > dmax) {
                dmax = d;
                closestFaceB = face;
              }
            }

            const worldVertsB1 = [];

            for (let i = 0; i < hullB.faces[closestFaceB].length; i++) {
              const b = hullB.vertices[hullB.faces[closestFaceB][i]];
              const worldb = new Vec3();
              worldb.copy(b);
              quatB.vmult(worldb, worldb);
              posB.vadd(worldb, worldb);
              worldVertsB1.push(worldb);
            }

            if (closestFaceB >= 0) {
              this.clipFaceAgainstHull(
                separatingNormal,
                posA,
                quatA,
                worldVertsB1,
                minDist,
                maxDist,
                result
              );
            }
          }
          /**
           * Find the separating axis between this hull and another
           * @param target The target vector to save the axis in
           * @return Returns false if a separation is found, else true
           */

          findSeparatingAxis(
            hullB,
            posA,
            quatA,
            posB,
            quatB,
            target,
            faceListA,
            faceListB
          ) {
            const faceANormalWS3 = new Vec3();
            const Worldnormal1 = new Vec3();
            const deltaC = new Vec3();
            const worldEdge0 = new Vec3();
            const worldEdge1 = new Vec3();
            const Cross = new Vec3();
            let dmin = Number.MAX_VALUE;
            const hullA = this;

            if (!hullA.uniqueAxes) {
              const numFacesA = faceListA
                ? faceListA.length
                : hullA.faces.length; // Test face normals from hullA

              for (let i = 0; i < numFacesA; i++) {
                const fi = faceListA ? faceListA[i] : i; // Get world face normal

                faceANormalWS3.copy(hullA.faceNormals[fi]);
                quatA.vmult(faceANormalWS3, faceANormalWS3);
                const d = hullA.testSepAxis(
                  faceANormalWS3,
                  hullB,
                  posA,
                  quatA,
                  posB,
                  quatB
                );

                if (d === false) {
                  return false;
                }

                if (d < dmin) {
                  dmin = d;
                  target.copy(faceANormalWS3);
                }
              }
            } else {
              // Test unique axes
              for (let i = 0; i !== hullA.uniqueAxes.length; i++) {
                // Get world axis
                quatA.vmult(hullA.uniqueAxes[i], faceANormalWS3);
                const d = hullA.testSepAxis(
                  faceANormalWS3,
                  hullB,
                  posA,
                  quatA,
                  posB,
                  quatB
                );

                if (d === false) {
                  return false;
                }

                if (d < dmin) {
                  dmin = d;
                  target.copy(faceANormalWS3);
                }
              }
            }

            if (!hullB.uniqueAxes) {
              // Test face normals from hullB
              const numFacesB = faceListB
                ? faceListB.length
                : hullB.faces.length;

              for (let i = 0; i < numFacesB; i++) {
                const fi = faceListB ? faceListB[i] : i;
                Worldnormal1.copy(hullB.faceNormals[fi]);
                quatB.vmult(Worldnormal1, Worldnormal1);
                const d = hullA.testSepAxis(
                  Worldnormal1,
                  hullB,
                  posA,
                  quatA,
                  posB,
                  quatB
                );

                if (d === false) {
                  return false;
                }

                if (d < dmin) {
                  dmin = d;
                  target.copy(Worldnormal1);
                }
              }
            } else {
              // Test unique axes in B
              for (let i = 0; i !== hullB.uniqueAxes.length; i++) {
                quatB.vmult(hullB.uniqueAxes[i], Worldnormal1);
                const d = hullA.testSepAxis(
                  Worldnormal1,
                  hullB,
                  posA,
                  quatA,
                  posB,
                  quatB
                );

                if (d === false) {
                  return false;
                }

                if (d < dmin) {
                  dmin = d;
                  target.copy(Worldnormal1);
                }
              }
            } // Test edges

            for (let e0 = 0; e0 !== hullA.uniqueEdges.length; e0++) {
              // Get world edge
              quatA.vmult(hullA.uniqueEdges[e0], worldEdge0);

              for (let e1 = 0; e1 !== hullB.uniqueEdges.length; e1++) {
                // Get world edge 2
                quatB.vmult(hullB.uniqueEdges[e1], worldEdge1);
                worldEdge0.cross(worldEdge1, Cross);

                if (!Cross.almostZero()) {
                  Cross.normalize();
                  const dist = hullA.testSepAxis(
                    Cross,
                    hullB,
                    posA,
                    quatA,
                    posB,
                    quatB
                  );

                  if (dist === false) {
                    return false;
                  }

                  if (dist < dmin) {
                    dmin = dist;
                    target.copy(Cross);
                  }
                }
              }
            }

            posB.vsub(posA, deltaC);

            if (deltaC.dot(target) > 0.0) {
              target.negate(target);
            }

            return true;
          }
          /**
           * Test separating axis against two hulls. Both hulls are projected onto the axis and the overlap size is returned if there is one.
           * @return The overlap depth, or FALSE if no penetration.
           */

          testSepAxis(axis, hullB, posA, quatA, posB, quatB) {
            const hullA = this;
            ConvexPolyhedron.project(hullA, axis, posA, quatA, maxminA);
            ConvexPolyhedron.project(hullB, axis, posB, quatB, maxminB);
            const maxA = maxminA[0];
            const minA = maxminA[1];
            const maxB = maxminB[0];
            const minB = maxminB[1];

            if (maxA < minB || maxB < minA) {
              return false; // Separated
            }

            const d0 = maxA - minB;
            const d1 = maxB - minA;
            const depth = d0 < d1 ? d0 : d1;
            return depth;
          }
          /**
           * calculateLocalInertia
           */

          calculateLocalInertia(mass, target) {
            // Approximate with box inertia
            // Exact inertia calculation is overkill, but see http://geometrictools.com/Documentation/PolyhedralMassProperties.pdf for the correct way to do it
            const aabbmax = new Vec3();
            const aabbmin = new Vec3();
            this.computeLocalAABB(aabbmin, aabbmax);
            const x = aabbmax.x - aabbmin.x;
            const y = aabbmax.y - aabbmin.y;
            const z = aabbmax.z - aabbmin.z;
            target.x = (1.0 / 12.0) * mass * (2 * y * 2 * y + 2 * z * 2 * z);
            target.y = (1.0 / 12.0) * mass * (2 * x * 2 * x + 2 * z * 2 * z);
            target.z = (1.0 / 12.0) * mass * (2 * y * 2 * y + 2 * x * 2 * x);
          }
          /**
           * @param face_i Index of the face
           */

          getPlaneConstantOfFace(face_i) {
            const f = this.faces[face_i];
            const n = this.faceNormals[face_i];
            const v = this.vertices[f[0]];
            const c = -n.dot(v);
            return c;
          }
          /**
           * Clip a face against a hull.
           * @param worldVertsB1 An array of Vec3 with vertices in the world frame.
           * @param minDist Distance clamping
           * @param Array result Array to store resulting contact points in. Will be objects with properties: point, depth, normal. These are represented in world coordinates.
           */

          clipFaceAgainstHull(
            separatingNormal,
            posA,
            quatA,
            worldVertsB1,
            minDist,
            maxDist,
            result
          ) {
            const faceANormalWS = new Vec3();
            const edge0 = new Vec3();
            const WorldEdge0 = new Vec3();
            const worldPlaneAnormal1 = new Vec3();
            const planeNormalWS1 = new Vec3();
            const worldA1 = new Vec3();
            const localPlaneNormal = new Vec3();
            const planeNormalWS = new Vec3();
            const hullA = this;
            const worldVertsB2 = [];
            const pVtxIn = worldVertsB1;
            const pVtxOut = worldVertsB2;
            let closestFaceA = -1;
            let dmin = Number.MAX_VALUE; // Find the face with normal closest to the separating axis

            for (let face = 0; face < hullA.faces.length; face++) {
              faceANormalWS.copy(hullA.faceNormals[face]);
              quatA.vmult(faceANormalWS, faceANormalWS);
              const d = faceANormalWS.dot(separatingNormal);

              if (d < dmin) {
                dmin = d;
                closestFaceA = face;
              }
            }

            if (closestFaceA < 0) {
              return;
            } // Get the face and construct connected faces

            const polyA = hullA.faces[closestFaceA];
            polyA.connectedFaces = [];

            for (let i = 0; i < hullA.faces.length; i++) {
              for (let j = 0; j < hullA.faces[i].length; j++) {
                if (
                  /* Sharing a vertex*/
                  polyA.indexOf(hullA.faces[i][j]) !== -1 &&
                  /* Not the one we are looking for connections from */
                  i !== closestFaceA &&
                  /* Not already added */
                  polyA.connectedFaces.indexOf(i) === -1
                ) {
                  polyA.connectedFaces.push(i);
                }
              }
            } // Clip the polygon to the back of the planes of all faces of hull A,
            // that are adjacent to the witness face

            const numVerticesA = polyA.length;

            for (let i = 0; i < numVerticesA; i++) {
              const a = hullA.vertices[polyA[i]];
              const b = hullA.vertices[polyA[(i + 1) % numVerticesA]];
              a.vsub(b, edge0);
              WorldEdge0.copy(edge0);
              quatA.vmult(WorldEdge0, WorldEdge0);
              posA.vadd(WorldEdge0, WorldEdge0);
              worldPlaneAnormal1.copy(this.faceNormals[closestFaceA]);
              quatA.vmult(worldPlaneAnormal1, worldPlaneAnormal1);
              posA.vadd(worldPlaneAnormal1, worldPlaneAnormal1);
              WorldEdge0.cross(worldPlaneAnormal1, planeNormalWS1);
              planeNormalWS1.negate(planeNormalWS1);
              worldA1.copy(a);
              quatA.vmult(worldA1, worldA1);
              posA.vadd(worldA1, worldA1);
              const otherFace = polyA.connectedFaces[i];
              localPlaneNormal.copy(this.faceNormals[otherFace]);
              const localPlaneEq = this.getPlaneConstantOfFace(otherFace);
              planeNormalWS.copy(localPlaneNormal);
              quatA.vmult(planeNormalWS, planeNormalWS);
              const planeEqWS = localPlaneEq - planeNormalWS.dot(posA); // Clip face against our constructed plane

              this.clipFaceAgainstPlane(
                pVtxIn,
                pVtxOut,
                planeNormalWS,
                planeEqWS
              ); // Throw away all clipped points, but save the remaining until next clip

              while (pVtxIn.length) {
                pVtxIn.shift();
              }

              while (pVtxOut.length) {
                pVtxIn.push(pVtxOut.shift());
              }
            } // only keep contact points that are behind the witness face

            localPlaneNormal.copy(this.faceNormals[closestFaceA]);
            const localPlaneEq = this.getPlaneConstantOfFace(closestFaceA);
            planeNormalWS.copy(localPlaneNormal);
            quatA.vmult(planeNormalWS, planeNormalWS);
            const planeEqWS = localPlaneEq - planeNormalWS.dot(posA);

            for (let i = 0; i < pVtxIn.length; i++) {
              let depth = planeNormalWS.dot(pVtxIn[i]) + planeEqWS; // ???

              if (depth <= minDist) {
                console.log(
                  "clamped: depth=" + depth + " to minDist=" + minDist
                );
                depth = minDist;
              }

              if (depth <= maxDist) {
                const point = pVtxIn[i];

                if (depth <= 1e-6) {
                  const p = {
                    point,
                    normal: planeNormalWS,
                    depth,
                  };
                  result.push(p);
                }
              }
            }
          }
          /**
           * Clip a face in a hull against the back of a plane.
           * @param planeConstant The constant in the mathematical plane equation
           */

          clipFaceAgainstPlane(
            inVertices,
            outVertices,
            planeNormal,
            planeConstant
          ) {
            let n_dot_first;
            let n_dot_last;
            const numVerts = inVertices.length;

            if (numVerts < 2) {
              return outVertices;
            }

            let firstVertex = inVertices[inVertices.length - 1];
            let lastVertex = inVertices[0];
            n_dot_first = planeNormal.dot(firstVertex) + planeConstant;

            for (let vi = 0; vi < numVerts; vi++) {
              lastVertex = inVertices[vi];
              n_dot_last = planeNormal.dot(lastVertex) + planeConstant;

              if (n_dot_first < 0) {
                if (n_dot_last < 0) {
                  // Start < 0, end < 0, so output lastVertex
                  const newv = new Vec3();
                  newv.copy(lastVertex);
                  outVertices.push(newv);
                } else {
                  // Start < 0, end >= 0, so output intersection
                  const newv = new Vec3();
                  firstVertex.lerp(
                    lastVertex,
                    n_dot_first / (n_dot_first - n_dot_last),
                    newv
                  );
                  outVertices.push(newv);
                }
              } else {
                if (n_dot_last < 0) {
                  // Start >= 0, end < 0 so output intersection and end
                  const newv = new Vec3();
                  firstVertex.lerp(
                    lastVertex,
                    n_dot_first / (n_dot_first - n_dot_last),
                    newv
                  );
                  outVertices.push(newv);
                  outVertices.push(lastVertex);
                }
              }

              firstVertex = lastVertex;
              n_dot_first = n_dot_last;
            }

            return outVertices;
          }
          /**
           * Updates `.worldVertices` and sets `.worldVerticesNeedsUpdate` to false.
           */

          computeWorldVertices(position, quat) {
            while (this.worldVertices.length < this.vertices.length) {
              this.worldVertices.push(new Vec3());
            }

            const verts = this.vertices;
            const worldVerts = this.worldVertices;

            for (let i = 0; i !== this.vertices.length; i++) {
              quat.vmult(verts[i], worldVerts[i]);
              position.vadd(worldVerts[i], worldVerts[i]);
            }

            this.worldVerticesNeedsUpdate = false;
          }

          computeLocalAABB(aabbmin, aabbmax) {
            const vertices = this.vertices;
            aabbmin.set(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE);
            aabbmax.set(
              -Number.MAX_VALUE,
              -Number.MAX_VALUE,
              -Number.MAX_VALUE
            );

            for (let i = 0; i < this.vertices.length; i++) {
              const v = vertices[i];

              if (v.x < aabbmin.x) {
                aabbmin.x = v.x;
              } else if (v.x > aabbmax.x) {
                aabbmax.x = v.x;
              }

              if (v.y < aabbmin.y) {
                aabbmin.y = v.y;
              } else if (v.y > aabbmax.y) {
                aabbmax.y = v.y;
              }

              if (v.z < aabbmin.z) {
                aabbmin.z = v.z;
              } else if (v.z > aabbmax.z) {
                aabbmax.z = v.z;
              }
            }
          }
          /**
           * Updates `worldVertices` and sets `worldVerticesNeedsUpdate` to false.
           */

          computeWorldFaceNormals(quat) {
            const N = this.faceNormals.length;

            while (this.worldFaceNormals.length < N) {
              this.worldFaceNormals.push(new Vec3());
            }

            const normals = this.faceNormals;
            const worldNormals = this.worldFaceNormals;

            for (let i = 0; i !== N; i++) {
              quat.vmult(normals[i], worldNormals[i]);
            }

            this.worldFaceNormalsNeedsUpdate = false;
          }
          /**
           * updateBoundingSphereRadius
           */

          updateBoundingSphereRadius() {
            // Assume points are distributed with local (0,0,0) as center
            let max2 = 0;
            const verts = this.vertices;

            for (let i = 0; i !== verts.length; i++) {
              const norm2 = verts[i].lengthSquared();

              if (norm2 > max2) {
                max2 = norm2;
              }
            }

            this.boundingSphereRadius = Math.sqrt(max2);
          }
          /**
           * calculateWorldAABB
           */

          calculateWorldAABB(pos, quat, min, max) {
            const verts = this.vertices;
            let minx;
            let miny;
            let minz;
            let maxx;
            let maxy;
            let maxz;
            let tempWorldVertex = new Vec3();

            for (let i = 0; i < verts.length; i++) {
              tempWorldVertex.copy(verts[i]);
              quat.vmult(tempWorldVertex, tempWorldVertex);
              pos.vadd(tempWorldVertex, tempWorldVertex);
              const v = tempWorldVertex;

              if (minx === undefined || v.x < minx) {
                minx = v.x;
              }

              if (maxx === undefined || v.x > maxx) {
                maxx = v.x;
              }

              if (miny === undefined || v.y < miny) {
                miny = v.y;
              }

              if (maxy === undefined || v.y > maxy) {
                maxy = v.y;
              }

              if (minz === undefined || v.z < minz) {
                minz = v.z;
              }

              if (maxz === undefined || v.z > maxz) {
                maxz = v.z;
              }
            }

            min.set(minx, miny, minz);
            max.set(maxx, maxy, maxz);
          }
          /**
           * Get approximate convex volume
           */

          volume() {
            return (4.0 * Math.PI * this.boundingSphereRadius) / 3.0;
          }
          /**
           * Get an average of all the vertices positions
           */

          getAveragePointLocal(target = new Vec3()) {
            const verts = this.vertices;

            for (let i = 0; i < verts.length; i++) {
              target.vadd(verts[i], target);
            }

            target.scale(1 / verts.length, target);
            return target;
          }
          /**
           * Transform all local points. Will change the .vertices
           */

          transformAllPoints(offset, quat) {
            const n = this.vertices.length;
            const verts = this.vertices; // Apply rotation

            if (quat) {
              // Rotate vertices
              for (let i = 0; i < n; i++) {
                const v = verts[i];
                quat.vmult(v, v);
              } // Rotate face normals

              for (let i = 0; i < this.faceNormals.length; i++) {
                const v = this.faceNormals[i];
                quat.vmult(v, v);
              }
              /*
            // Rotate edges
            for(let i=0; i<this.uniqueEdges.length; i++){
                const v = this.uniqueEdges[i];
                quat.vmult(v,v);
            }*/
            } // Apply offset

            if (offset) {
              for (let i = 0; i < n; i++) {
                const v = verts[i];
                v.vadd(offset, v);
              }
            }
          }
          /**
           * Checks whether p is inside the polyhedra. Must be in local coords.
           * The point lies outside of the convex hull of the other points if and only if the direction
           * of all the vectors from it to those other points are on less than one half of a sphere around it.
           * @param p A point given in local coordinates
           */

          pointIsInside(p) {
            const verts = this.vertices;
            const faces = this.faces;
            const normals = this.faceNormals;
            const pointInside = new Vec3();
            this.getAveragePointLocal(pointInside);

            for (let i = 0; i < this.faces.length; i++) {
              let n = normals[i];
              const v = verts[faces[i][0]]; // We only need one point in the face
              // This dot product determines which side of the edge the point is

              const vToP = new Vec3();
              p.vsub(v, vToP);
              const r1 = n.dot(vToP);
              const vToPointInside = new Vec3();
              pointInside.vsub(v, vToPointInside);
              const r2 = n.dot(vToPointInside);

              if ((r1 < 0 && r2 > 0) || (r1 > 0 && r2 < 0)) {
                return false; // Encountered some other sign. Exit.
              }
            } // If we got here, all dot products were of the same sign.

            return -1;
          }
          /**
           * Get max and min dot product of a convex hull at position (pos,quat) projected onto an axis.
           * Results are saved in the array maxmin.
           * @param result result[0] and result[1] will be set to maximum and minimum, respectively.
           */

          static project(shape, axis, pos, quat, result) {
            const n = shape.vertices.length;
            const localAxis = project_localAxis;
            let max = 0;
            let min = 0;
            const localOrigin = project_localOrigin;
            const vs = shape.vertices;
            localOrigin.setZero(); // Transform the axis to local

            Transform.vectorToLocalFrame(pos, quat, axis, localAxis);
            Transform.pointToLocalFrame(pos, quat, localOrigin, localOrigin);
            const add = localOrigin.dot(localAxis);
            min = max = vs[0].dot(localAxis);

            for (let i = 1; i < n; i++) {
              const val = vs[i].dot(localAxis);

              if (val > max) {
                max = val;
              }

              if (val < min) {
                min = val;
              }
            }

            min -= add;
            max -= add;

            if (min > max) {
              // Inconsistent - swap
              const temp = min;
              min = max;
              max = temp;
            } // Output

            result[0] = max;
            result[1] = min;
          }
        }
        const maxminA = [];
        const maxminB = [];
        const project_localAxis = new Vec3();
        const project_localOrigin = new Vec3();

        /**
         * A 3d box shape.
         * @example
         *     const size = 1
         *     const halfExtents = new CANNON.Vec3(size, size, size)
         *     const boxShape = new CANNON.Box(halfExtents)
         *     const boxBody = new CANNON.Body({ mass: 1, shape: boxShape })
         *     world.addBody(boxBody)
         */
        class Box extends Shape {
          /**
           * The half extents of the box.
           */

          /**
           * Used by the contact generator to make contacts with other convex polyhedra for example.
           */
          constructor(halfExtents) {
            super({
              type: Shape.types.BOX,
            });
            this.halfExtents = void 0;
            this.convexPolyhedronRepresentation = void 0;
            this.halfExtents = halfExtents;
            this.convexPolyhedronRepresentation = null;
            this.updateConvexPolyhedronRepresentation();
            this.updateBoundingSphereRadius();
          }
          /**
           * Updates the local convex polyhedron representation used for some collisions.
           */

          updateConvexPolyhedronRepresentation() {
            const sx = this.halfExtents.x;
            const sy = this.halfExtents.y;
            const sz = this.halfExtents.z;
            const V = Vec3;
            const vertices = [
              new V(-sx, -sy, -sz),
              new V(sx, -sy, -sz),
              new V(sx, sy, -sz),
              new V(-sx, sy, -sz),
              new V(-sx, -sy, sz),
              new V(sx, -sy, sz),
              new V(sx, sy, sz),
              new V(-sx, sy, sz),
            ];
            const faces = [
              [3, 2, 1, 0], // -z
              [4, 5, 6, 7], // +z
              [5, 4, 0, 1], // -y
              [2, 3, 7, 6], // +y
              [0, 4, 7, 3], // -x
              [1, 2, 6, 5], // +x
            ];
            const axes = [new V(0, 0, 1), new V(0, 1, 0), new V(1, 0, 0)];
            const h = new ConvexPolyhedron({
              vertices,
              faces,
              axes,
            });
            this.convexPolyhedronRepresentation = h;
            h.material = this.material;
          }
          /**
           * Calculate the inertia of the box.
           */

          calculateLocalInertia(mass, target = new Vec3()) {
            Box.calculateInertia(this.halfExtents, mass, target);
            return target;
          }

          static calculateInertia(halfExtents, mass, target) {
            const e = halfExtents;
            target.x =
              (1.0 / 12.0) * mass * (2 * e.y * 2 * e.y + 2 * e.z * 2 * e.z);
            target.y =
              (1.0 / 12.0) * mass * (2 * e.x * 2 * e.x + 2 * e.z * 2 * e.z);
            target.z =
              (1.0 / 12.0) * mass * (2 * e.y * 2 * e.y + 2 * e.x * 2 * e.x);
          }
          /**
           * Get the box 6 side normals
           * @param sixTargetVectors An array of 6 vectors, to store the resulting side normals in.
           * @param quat Orientation to apply to the normal vectors. If not provided, the vectors will be in respect to the local frame.
           */

          getSideNormals(sixTargetVectors, quat) {
            const sides = sixTargetVectors;
            const ex = this.halfExtents;
            sides[0].set(ex.x, 0, 0);
            sides[1].set(0, ex.y, 0);
            sides[2].set(0, 0, ex.z);
            sides[3].set(-ex.x, 0, 0);
            sides[4].set(0, -ex.y, 0);
            sides[5].set(0, 0, -ex.z);

            if (quat !== undefined) {
              for (let i = 0; i !== sides.length; i++) {
                quat.vmult(sides[i], sides[i]);
              }
            }

            return sides;
          }
          /**
           * Returns the volume of the box.
           */

          volume() {
            return (
              8.0 * this.halfExtents.x * this.halfExtents.y * this.halfExtents.z
            );
          }
          /**
           * updateBoundingSphereRadius
           */

          updateBoundingSphereRadius() {
            this.boundingSphereRadius = this.halfExtents.length();
          }
          /**
           * forEachWorldCorner
           */

          forEachWorldCorner(pos, quat, callback) {
            const e = this.halfExtents;
            const corners = [
              [e.x, e.y, e.z],
              [-e.x, e.y, e.z],
              [-e.x, -e.y, e.z],
              [-e.x, -e.y, -e.z],
              [e.x, -e.y, -e.z],
              [e.x, e.y, -e.z],
              [-e.x, e.y, -e.z],
              [e.x, -e.y, e.z],
            ];

            for (let i = 0; i < corners.length; i++) {
              worldCornerTempPos.set(
                corners[i][0],
                corners[i][1],
                corners[i][2]
              );
              quat.vmult(worldCornerTempPos, worldCornerTempPos);
              pos.vadd(worldCornerTempPos, worldCornerTempPos);
              callback(
                worldCornerTempPos.x,
                worldCornerTempPos.y,
                worldCornerTempPos.z
              );
            }
          }
          /**
           * calculateWorldAABB
           */

          calculateWorldAABB(pos, quat, min, max) {
            const e = this.halfExtents;
            worldCornersTemp[0].set(e.x, e.y, e.z);
            worldCornersTemp[1].set(-e.x, e.y, e.z);
            worldCornersTemp[2].set(-e.x, -e.y, e.z);
            worldCornersTemp[3].set(-e.x, -e.y, -e.z);
            worldCornersTemp[4].set(e.x, -e.y, -e.z);
            worldCornersTemp[5].set(e.x, e.y, -e.z);
            worldCornersTemp[6].set(-e.x, e.y, -e.z);
            worldCornersTemp[7].set(e.x, -e.y, e.z);
            const wc = worldCornersTemp[0];
            quat.vmult(wc, wc);
            pos.vadd(wc, wc);
            max.copy(wc);
            min.copy(wc);

            for (let i = 1; i < 8; i++) {
              const wc = worldCornersTemp[i];
              quat.vmult(wc, wc);
              pos.vadd(wc, wc);
              const x = wc.x;
              const y = wc.y;
              const z = wc.z;

              if (x > max.x) {
                max.x = x;
              }

              if (y > max.y) {
                max.y = y;
              }

              if (z > max.z) {
                max.z = z;
              }

              if (x < min.x) {
                min.x = x;
              }

              if (y < min.y) {
                min.y = y;
              }

              if (z < min.z) {
                min.z = z;
              }
            } // Get each axis max
            // min.set(Infinity,Infinity,Infinity);
            // max.set(-Infinity,-Infinity,-Infinity);
            // this.forEachWorldCorner(pos,quat,function(x,y,z){
            //     if(x > max.x){
            //         max.x = x;
            //     }
            //     if(y > max.y){
            //         max.y = y;
            //     }
            //     if(z > max.z){
            //         max.z = z;
            //     }
            //     if(x < min.x){
            //         min.x = x;
            //     }
            //     if(y < min.y){
            //         min.y = y;
            //     }
            //     if(z < min.z){
            //         min.z = z;
            //     }
            // });
          }
        }
        const worldCornerTempPos = new Vec3();
        const worldCornersTemp = [
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
          new Vec3(),
        ];

        /**
         * BODY_TYPES
         */
        const BODY_TYPES = {
          /** DYNAMIC */
          DYNAMIC: 1,

          /** STATIC */
          STATIC: 2,

          /** KINEMATIC */
          KINEMATIC: 4,
        };
        /**
         * BodyType
         */

        /**
         * BODY_SLEEP_STATES
         */
        const BODY_SLEEP_STATES = {
          /** AWAKE */
          AWAKE: 0,

          /** SLEEPY */
          SLEEPY: 1,

          /** SLEEPING */
          SLEEPING: 2,
        };
        /**
         * BodySleepState
         */

        /**
         * Base class for all body types.
         * @example
         *     const shape = new CANNON.Sphere(1)
         *     const body = new CANNON.Body({
         *       mass: 1,
         *       shape,
         *     })
         *     world.addBody(body)
         */
        class Body extends EventTarget {
          /**
           * Dispatched after two bodies collide. This event is dispatched on each
           * of the two bodies involved in the collision.
           * @event collide
           * @param body The body that was involved in the collision.
           * @param contact The details of the collision.
           */

          /**
           * A dynamic body is fully simulated. Can be moved manually by the user, but normally they move according to forces. A dynamic body can collide with all body types. A dynamic body always has finite, non-zero mass.
           */

          /**
           * A static body does not move during simulation and behaves as if it has infinite mass. Static bodies can be moved manually by setting the position of the body. The velocity of a static body is always zero. Static bodies do not collide with other static or kinematic bodies.
           */

          /**
           * A kinematic body moves under simulation according to its velocity. They do not respond to forces. They can be moved manually, but normally a kinematic body is moved by setting its velocity. A kinematic body behaves as if it has infinite mass. Kinematic bodies do not collide with other static or kinematic bodies.
           */

          /**
           * AWAKE
           */

          /**
           * SLEEPY
           */

          /**
           * SLEEPING
           */

          /**
           * Dispatched after a sleeping body has woken up.
           * @event wakeup
           */

          /**
           * Dispatched after a body has gone in to the sleepy state.
           * @event sleepy
           */

          /**
           * Dispatched after a body has fallen asleep.
           * @event sleep
           */

          /**
           * Identifier of the body.
           */

          /**
           * Position of body in World.bodies. Updated by World and used in ArrayCollisionMatrix.
           */

          /**
           * Reference to the world the body is living in.
           */

          /**
           * Callback function that is used BEFORE stepping the system. Use it to apply forces, for example. Inside the function, "this" will refer to this Body object. Deprecated - use World events instead.
           * @deprecated Use World events instead
           */

          /**
           * Callback function that is used AFTER stepping the system. Inside the function, "this" will refer to this Body object. Deprecated - use World events instead.
           * @deprecated Use World events instead
           */

          /**
           * The collision group the body belongs to.
           * @default 1
           */

          /**
           * The collision group the body can collide with.
           * @default -1
           */

          /**
           * Whether to produce contact forces when in contact with other bodies. Note that contacts will be generated, but they will be disabled - i.e. "collide" events will be raised, but forces will not be altered.
           */

          /**
           * World space position of the body.
           */

          /**
           * Interpolated position of the body.
           */

          /**
           * Initial position of the body.
           */

          /**
           * World space velocity of the body.
           */

          /**
           * Initial velocity of the body.
           */

          /**
           * Linear force on the body in world space.
           */

          /**
           * The mass of the body.
           * @default 0
           */

          /**
           * The physics material of the body. It defines the body interaction with other bodies.
           */

          /**
           * How much to damp the body velocity each step. It can go from 0 to 1.
           * @default 0.01
           */

          /**
           * One of: `Body.DYNAMIC`, `Body.STATIC` and `Body.KINEMATIC`.
           */

          /**
           * If true, the body will automatically fall to sleep.
           * @default true
           */

          /**
           * Current sleep state.
           */

          /**
           * If the speed (the norm of the velocity) is smaller than this value, the body is considered sleepy.
           * @default 0.1
           */

          /**
           * If the body has been sleepy for this sleepTimeLimit seconds, it is considered sleeping.
           * @default 1
           */

          /**
           * World space rotational force on the body, around center of mass.
           */

          /**
           * World space orientation of the body.
           */

          /**
           * Initial quaternion of the body.
           */

          /**
           * Interpolated orientation of the body.
           */

          /**
           * Angular velocity of the body, in world space. Think of the angular velocity as a vector, which the body rotates around. The length of this vector determines how fast (in radians per second) the body rotates.
           */

          /**
           * Initial angular velocity of the body.
           */

          /**
           * List of Shapes that have been added to the body.
           */

          /**
           * Position of each Shape in the body, given in local Body space.
           */

          /**
           * Orientation of each Shape, given in local Body space.
           */

          /**
           * The inertia of the body.
           */

          /**
           * Set to true if you don't want the body to rotate. Make sure to run .updateMassProperties() if you change this after the body creation.
           * @default false
           */

          /**
           * How much to damp the body angular velocity each step. It can go from 0 to 1.
           * @default 0.01
           */

          /**
           * Use this property to limit the motion along any world axis. (1,1,1) will allow motion along all axes while (0,0,0) allows none.
           */

          /**
           * Use this property to limit the rotational motion along any world axis. (1,1,1) will allow rotation along all axes while (0,0,0) allows none.
           */

          /**
           * World space bounding box of the body and its shapes.
           */

          /**
           * Indicates if the AABB needs to be updated before use.
           */

          /**
           * Total bounding radius of the Body including its shapes, relative to body.position.
           */

          /**
           * When true the body behaves like a trigger. It does not collide
           * with other bodies but collision events are still triggered.
           * @default false
           */
          constructor(options = {}) {
            super();
            this.id = void 0;
            this.index = void 0;
            this.world = void 0;
            this.preStep = void 0;
            this.postStep = void 0;
            this.vlambda = void 0;
            this.collisionFilterGroup = void 0;
            this.collisionFilterMask = void 0;
            this.collisionResponse = void 0;
            this.position = void 0;
            this.previousPosition = void 0;
            this.interpolatedPosition = void 0;
            this.initPosition = void 0;
            this.velocity = void 0;
            this.initVelocity = void 0;
            this.force = void 0;
            this.mass = void 0;
            this.invMass = void 0;
            this.material = void 0;
            this.linearDamping = void 0;
            this.type = void 0;
            this.allowSleep = void 0;
            this.sleepState = void 0;
            this.sleepSpeedLimit = void 0;
            this.sleepTimeLimit = void 0;
            this.timeLastSleepy = void 0;
            this.wakeUpAfterNarrowphase = void 0;
            this.torque = void 0;
            this.quaternion = void 0;
            this.initQuaternion = void 0;
            this.previousQuaternion = void 0;
            this.interpolatedQuaternion = void 0;
            this.angularVelocity = void 0;
            this.initAngularVelocity = void 0;
            this.shapes = void 0;
            this.shapeOffsets = void 0;
            this.shapeOrientations = void 0;
            this.inertia = void 0;
            this.invInertia = void 0;
            this.invInertiaWorld = void 0;
            this.invMassSolve = void 0;
            this.invInertiaSolve = void 0;
            this.invInertiaWorldSolve = void 0;
            this.fixedRotation = void 0;
            this.angularDamping = void 0;
            this.linearFactor = void 0;
            this.angularFactor = void 0;
            this.aabb = void 0;
            this.aabbNeedsUpdate = void 0;
            this.boundingRadius = void 0;
            this.wlambda = void 0;
            this.isTrigger = void 0;
            this.id = Body.idCounter++;
            this.index = -1;
            this.world = null;
            this.preStep = null;
            this.postStep = null;
            this.vlambda = new Vec3();
            this.collisionFilterGroup =
              typeof options.collisionFilterGroup === "number"
                ? options.collisionFilterGroup
                : 1;
            this.collisionFilterMask =
              typeof options.collisionFilterMask === "number"
                ? options.collisionFilterMask
                : -1;
            this.collisionResponse =
              typeof options.collisionResponse === "boolean"
                ? options.collisionResponse
                : true;
            this.position = new Vec3();
            this.previousPosition = new Vec3();
            this.interpolatedPosition = new Vec3();
            this.initPosition = new Vec3();

            if (options.position) {
              this.position.copy(options.position);
              this.previousPosition.copy(options.position);
              this.interpolatedPosition.copy(options.position);
              this.initPosition.copy(options.position);
            }

            this.velocity = new Vec3();

            if (options.velocity) {
              this.velocity.copy(options.velocity);
            }

            this.initVelocity = new Vec3();
            this.force = new Vec3();
            const mass = typeof options.mass === "number" ? options.mass : 0;
            this.mass = mass;
            this.invMass = mass > 0 ? 1.0 / mass : 0;
            this.material = options.material || null;
            this.linearDamping =
              typeof options.linearDamping === "number"
                ? options.linearDamping
                : 0.01;
            this.type = mass <= 0.0 ? Body.STATIC : Body.DYNAMIC;

            if (typeof options.type === typeof Body.STATIC) {
              this.type = options.type;
            }

            this.allowSleep =
              typeof options.allowSleep !== "undefined"
                ? options.allowSleep
                : true;
            this.sleepState = Body.AWAKE;
            this.sleepSpeedLimit =
              typeof options.sleepSpeedLimit !== "undefined"
                ? options.sleepSpeedLimit
                : 0.1;
            this.sleepTimeLimit =
              typeof options.sleepTimeLimit !== "undefined"
                ? options.sleepTimeLimit
                : 1;
            this.timeLastSleepy = 0;
            this.wakeUpAfterNarrowphase = false;
            this.torque = new Vec3();
            this.quaternion = new Quaternion();
            this.initQuaternion = new Quaternion();
            this.previousQuaternion = new Quaternion();
            this.interpolatedQuaternion = new Quaternion();

            if (options.quaternion) {
              this.quaternion.copy(options.quaternion);
              this.initQuaternion.copy(options.quaternion);
              this.previousQuaternion.copy(options.quaternion);
              this.interpolatedQuaternion.copy(options.quaternion);
            }

            this.angularVelocity = new Vec3();

            if (options.angularVelocity) {
              this.angularVelocity.copy(options.angularVelocity);
            }

            this.initAngularVelocity = new Vec3();
            this.shapes = [];
            this.shapeOffsets = [];
            this.shapeOrientations = [];
            this.inertia = new Vec3();
            this.invInertia = new Vec3();
            this.invInertiaWorld = new Mat3();
            this.invMassSolve = 0;
            this.invInertiaSolve = new Vec3();
            this.invInertiaWorldSolve = new Mat3();
            this.fixedRotation =
              typeof options.fixedRotation !== "undefined"
                ? options.fixedRotation
                : false;
            this.angularDamping =
              typeof options.angularDamping !== "undefined"
                ? options.angularDamping
                : 0.01;
            this.linearFactor = new Vec3(1, 1, 1);

            if (options.linearFactor) {
              this.linearFactor.copy(options.linearFactor);
            }

            this.angularFactor = new Vec3(1, 1, 1);

            if (options.angularFactor) {
              this.angularFactor.copy(options.angularFactor);
            }

            this.aabb = new AABB();
            this.aabbNeedsUpdate = true;
            this.boundingRadius = 0;
            this.wlambda = new Vec3();
            this.isTrigger = Boolean(options.isTrigger);

            if (options.shape) {
              this.addShape(options.shape);
            }

            this.updateMassProperties();
          }
          /**
           * Wake the body up.
           */

          wakeUp() {
            const prevState = this.sleepState;
            this.sleepState = Body.AWAKE;
            this.wakeUpAfterNarrowphase = false;

            if (prevState === Body.SLEEPING) {
              this.dispatchEvent(Body.wakeupEvent);
            }
          }
          /**
           * Force body sleep
           */

          sleep() {
            this.sleepState = Body.SLEEPING;
            this.velocity.set(0, 0, 0);
            this.angularVelocity.set(0, 0, 0);
            this.wakeUpAfterNarrowphase = false;
          }
          /**
           * Called every timestep to update internal sleep timer and change sleep state if needed.
           * @param time The world time in seconds
           */

          sleepTick(time) {
            if (this.allowSleep) {
              const sleepState = this.sleepState;
              const speedSquared =
                this.velocity.lengthSquared() +
                this.angularVelocity.lengthSquared();
              const speedLimitSquared = this.sleepSpeedLimit ** 2;

              if (
                sleepState === Body.AWAKE &&
                speedSquared < speedLimitSquared
              ) {
                this.sleepState = Body.SLEEPY; // Sleepy

                this.timeLastSleepy = time;
                this.dispatchEvent(Body.sleepyEvent);
              } else if (
                sleepState === Body.SLEEPY &&
                speedSquared > speedLimitSquared
              ) {
                this.wakeUp(); // Wake up
              } else if (
                sleepState === Body.SLEEPY &&
                time - this.timeLastSleepy > this.sleepTimeLimit
              ) {
                this.sleep(); // Sleeping

                this.dispatchEvent(Body.sleepEvent);
              }
            }
          }
          /**
           * If the body is sleeping, it should be immovable / have infinite mass during solve. We solve it by having a separate "solve mass".
           */

          updateSolveMassProperties() {
            if (
              this.sleepState === Body.SLEEPING ||
              this.type === Body.KINEMATIC
            ) {
              this.invMassSolve = 0;
              this.invInertiaSolve.setZero();
              this.invInertiaWorldSolve.setZero();
            } else {
              this.invMassSolve = this.invMass;
              this.invInertiaSolve.copy(this.invInertia);
              this.invInertiaWorldSolve.copy(this.invInertiaWorld);
            }
          }
          /**
           * Convert a world point to local body frame.
           */

          pointToLocalFrame(worldPoint, result = new Vec3()) {
            worldPoint.vsub(this.position, result);
            this.quaternion.conjugate().vmult(result, result);
            return result;
          }
          /**
           * Convert a world vector to local body frame.
           */

          vectorToLocalFrame(worldVector, result = new Vec3()) {
            this.quaternion.conjugate().vmult(worldVector, result);
            return result;
          }
          /**
           * Convert a local body point to world frame.
           */

          pointToWorldFrame(localPoint, result = new Vec3()) {
            this.quaternion.vmult(localPoint, result);
            result.vadd(this.position, result);
            return result;
          }
          /**
           * Convert a local body point to world frame.
           */

          vectorToWorldFrame(localVector, result = new Vec3()) {
            this.quaternion.vmult(localVector, result);
            return result;
          }
          /**
           * Add a shape to the body with a local offset and orientation.
           * @return The body object, for chainability.
           */

          addShape(shape, _offset, _orientation) {
            const offset = new Vec3();
            const orientation = new Quaternion();

            if (_offset) {
              offset.copy(_offset);
            }

            if (_orientation) {
              orientation.copy(_orientation);
            }

            this.shapes.push(shape);
            this.shapeOffsets.push(offset);
            this.shapeOrientations.push(orientation);
            this.updateMassProperties();
            this.updateBoundingRadius();
            this.aabbNeedsUpdate = true;
            shape.body = this;
            return this;
          }
          /**
           * Remove a shape from the body.
           * @return The body object, for chainability.
           */

          removeShape(shape) {
            const index = this.shapes.indexOf(shape);

            if (index === -1) {
              console.warn("Shape does not belong to the body");
              return this;
            }

            this.shapes.splice(index, 1);
            this.shapeOffsets.splice(index, 1);
            this.shapeOrientations.splice(index, 1);
            this.updateMassProperties();
            this.updateBoundingRadius();
            this.aabbNeedsUpdate = true;
            shape.body = null;
            return this;
          }
          /**
           * Update the bounding radius of the body. Should be done if any of the shapes are changed.
           */

          updateBoundingRadius() {
            const shapes = this.shapes;
            const shapeOffsets = this.shapeOffsets;
            const N = shapes.length;
            let radius = 0;

            for (let i = 0; i !== N; i++) {
              const shape = shapes[i];
              shape.updateBoundingSphereRadius();
              const offset = shapeOffsets[i].length();
              const r = shape.boundingSphereRadius;

              if (offset + r > radius) {
                radius = offset + r;
              }
            }

            this.boundingRadius = radius;
          }
          /**
           * Updates the .aabb
           */

          updateAABB() {
            const shapes = this.shapes;
            const shapeOffsets = this.shapeOffsets;
            const shapeOrientations = this.shapeOrientations;
            const N = shapes.length;
            const offset = tmpVec;
            const orientation = tmpQuat;
            const bodyQuat = this.quaternion;
            const aabb = this.aabb;
            const shapeAABB = updateAABB_shapeAABB;

            for (let i = 0; i !== N; i++) {
              const shape = shapes[i]; // Get shape world position

              bodyQuat.vmult(shapeOffsets[i], offset);
              offset.vadd(this.position, offset); // Get shape world quaternion

              bodyQuat.mult(shapeOrientations[i], orientation); // Get shape AABB

              shape.calculateWorldAABB(
                offset,
                orientation,
                shapeAABB.lowerBound,
                shapeAABB.upperBound
              );

              if (i === 0) {
                aabb.copy(shapeAABB);
              } else {
                aabb.extend(shapeAABB);
              }
            }

            this.aabbNeedsUpdate = false;
          }
          /**
           * Update `.inertiaWorld` and `.invInertiaWorld`
           */

          updateInertiaWorld(force) {
            const I = this.invInertia;

            if (I.x === I.y && I.y === I.z && !force);
            else {
              const m1 = uiw_m1;
              const m2 = uiw_m2;
              m1.setRotationFromQuaternion(this.quaternion);
              m1.transpose(m2);
              m1.scale(I, m1);
              m1.mmult(m2, this.invInertiaWorld);
            }
          }
          /**
           * Apply force to a point of the body. This could for example be a point on the Body surface.
           * Applying force this way will add to Body.force and Body.torque.
           * @param force The amount of force to add.
           * @param relativePoint A point relative to the center of mass to apply the force on.
           */

          applyForce(force, relativePoint = new Vec3()) {
            // Needed?
            if (this.type !== Body.DYNAMIC) {
              return;
            }

            if (this.sleepState === Body.SLEEPING) {
              this.wakeUp();
            } // Compute produced rotational force

            const rotForce = Body_applyForce_rotForce;
            relativePoint.cross(force, rotForce); // Add linear force

            this.force.vadd(force, this.force); // Add rotational force

            this.torque.vadd(rotForce, this.torque);
          }
          /**
           * Apply force to a local point in the body.
           * @param force The force vector to apply, defined locally in the body frame.
           * @param localPoint A local point in the body to apply the force on.
           */

          applyLocalForce(localForce, localPoint = new Vec3()) {
            if (this.type !== Body.DYNAMIC) {
              return;
            }

            const worldForce = Body_applyLocalForce_worldForce;
            const relativePointWorld = Body_applyLocalForce_relativePointWorld; // Transform the force vector to world space

            this.vectorToWorldFrame(localForce, worldForce);
            this.vectorToWorldFrame(localPoint, relativePointWorld);
            this.applyForce(worldForce, relativePointWorld);
          }
          /**
           * Apply torque to the body.
           * @param torque The amount of torque to add.
           */

          applyTorque(torque) {
            if (this.type !== Body.DYNAMIC) {
              return;
            }

            if (this.sleepState === Body.SLEEPING) {
              this.wakeUp();
            } // Add rotational force

            this.torque.vadd(torque, this.torque);
          }
          /**
           * Apply impulse to a point of the body. This could for example be a point on the Body surface.
           * An impulse is a force added to a body during a short period of time (impulse = force * time).
           * Impulses will be added to Body.velocity and Body.angularVelocity.
           * @param impulse The amount of impulse to add.
           * @param relativePoint A point relative to the center of mass to apply the force on.
           */

          applyImpulse(impulse, relativePoint = new Vec3()) {
            if (this.type !== Body.DYNAMIC) {
              return;
            }

            if (this.sleepState === Body.SLEEPING) {
              this.wakeUp();
            } // Compute point position relative to the body center

            const r = relativePoint; // Compute produced central impulse velocity

            const velo = Body_applyImpulse_velo;
            velo.copy(impulse);
            velo.scale(this.invMass, velo); // Add linear impulse

            this.velocity.vadd(velo, this.velocity); // Compute produced rotational impulse velocity

            const rotVelo = Body_applyImpulse_rotVelo;
            r.cross(impulse, rotVelo);
            /*
     rotVelo.x *= this.invInertia.x;
     rotVelo.y *= this.invInertia.y;
     rotVelo.z *= this.invInertia.z;
     */

            this.invInertiaWorld.vmult(rotVelo, rotVelo); // Add rotational Impulse

            this.angularVelocity.vadd(rotVelo, this.angularVelocity);
          }
          /**
           * Apply locally-defined impulse to a local point in the body.
           * @param force The force vector to apply, defined locally in the body frame.
           * @param localPoint A local point in the body to apply the force on.
           */

          applyLocalImpulse(localImpulse, localPoint = new Vec3()) {
            if (this.type !== Body.DYNAMIC) {
              return;
            }

            const worldImpulse = Body_applyLocalImpulse_worldImpulse;
            const relativePointWorld = Body_applyLocalImpulse_relativePoint; // Transform the force vector to world space

            this.vectorToWorldFrame(localImpulse, worldImpulse);
            this.vectorToWorldFrame(localPoint, relativePointWorld);
            this.applyImpulse(worldImpulse, relativePointWorld);
          }
          /**
           * Should be called whenever you change the body shape or mass.
           */

          updateMassProperties() {
            const halfExtents = Body_updateMassProperties_halfExtents;
            this.invMass = this.mass > 0 ? 1.0 / this.mass : 0;
            const I = this.inertia;
            const fixed = this.fixedRotation; // Approximate with AABB box

            this.updateAABB();
            halfExtents.set(
              (this.aabb.upperBound.x - this.aabb.lowerBound.x) / 2,
              (this.aabb.upperBound.y - this.aabb.lowerBound.y) / 2,
              (this.aabb.upperBound.z - this.aabb.lowerBound.z) / 2
            );
            Box.calculateInertia(halfExtents, this.mass, I);
            this.invInertia.set(
              I.x > 0 && !fixed ? 1.0 / I.x : 0,
              I.y > 0 && !fixed ? 1.0 / I.y : 0,
              I.z > 0 && !fixed ? 1.0 / I.z : 0
            );
            this.updateInertiaWorld(true);
          }
          /**
           * Get world velocity of a point in the body.
           * @param worldPoint
           * @param result
           * @return The result vector.
           */

          getVelocityAtWorldPoint(worldPoint, result) {
            const r = new Vec3();
            worldPoint.vsub(this.position, r);
            this.angularVelocity.cross(r, result);
            this.velocity.vadd(result, result);
            return result;
          }
          /**
           * Move the body forward in time.
           * @param dt Time step
           * @param quatNormalize Set to true to normalize the body quaternion
           * @param quatNormalizeFast If the quaternion should be normalized using "fast" quaternion normalization
           */

          integrate(dt, quatNormalize, quatNormalizeFast) {
            // Save previous position
            this.previousPosition.copy(this.position);
            this.previousQuaternion.copy(this.quaternion);

            if (
              !(this.type === Body.DYNAMIC || this.type === Body.KINEMATIC) ||
              this.sleepState === Body.SLEEPING
            ) {
              // Only for dynamic
              return;
            }

            const velo = this.velocity;
            const angularVelo = this.angularVelocity;
            const pos = this.position;
            const force = this.force;
            const torque = this.torque;
            const quat = this.quaternion;
            const invMass = this.invMass;
            const invInertia = this.invInertiaWorld;
            const linearFactor = this.linearFactor;
            const iMdt = invMass * dt;
            velo.x += force.x * iMdt * linearFactor.x;
            velo.y += force.y * iMdt * linearFactor.y;
            velo.z += force.z * iMdt * linearFactor.z;
            const e = invInertia.elements;
            const angularFactor = this.angularFactor;
            const tx = torque.x * angularFactor.x;
            const ty = torque.y * angularFactor.y;
            const tz = torque.z * angularFactor.z;
            angularVelo.x += dt * (e[0] * tx + e[1] * ty + e[2] * tz);
            angularVelo.y += dt * (e[3] * tx + e[4] * ty + e[5] * tz);
            angularVelo.z += dt * (e[6] * tx + e[7] * ty + e[8] * tz); // Use new velocity  - leap frog

            pos.x += velo.x * dt;
            pos.y += velo.y * dt;
            pos.z += velo.z * dt;
            quat.integrate(this.angularVelocity, dt, this.angularFactor, quat);

            if (quatNormalize) {
              if (quatNormalizeFast) {
                quat.normalizeFast();
              } else {
                quat.normalize();
              }
            }

            this.aabbNeedsUpdate = true; // Update world inertia

            this.updateInertiaWorld();
          }
        }
        Body.idCounter = 0;
        Body.COLLIDE_EVENT_NAME = "collide";
        Body.DYNAMIC = BODY_TYPES.DYNAMIC;
        Body.STATIC = BODY_TYPES.STATIC;
        Body.KINEMATIC = BODY_TYPES.KINEMATIC;
        Body.AWAKE = BODY_SLEEP_STATES.AWAKE;
        Body.SLEEPY = BODY_SLEEP_STATES.SLEEPY;
        Body.SLEEPING = BODY_SLEEP_STATES.SLEEPING;
        Body.wakeupEvent = {
          type: "wakeup",
        };
        Body.sleepyEvent = {
          type: "sleepy",
        };
        Body.sleepEvent = {
          type: "sleep",
        };
        const tmpVec = new Vec3();
        const tmpQuat = new Quaternion();
        const updateAABB_shapeAABB = new AABB();
        const uiw_m1 = new Mat3();
        const uiw_m2 = new Mat3();
        const Body_applyForce_rotForce = new Vec3();
        const Body_applyLocalForce_worldForce = new Vec3();
        const Body_applyLocalForce_relativePointWorld = new Vec3();
        const Body_applyImpulse_velo = new Vec3();
        const Body_applyImpulse_rotVelo = new Vec3();
        const Body_applyLocalImpulse_worldImpulse = new Vec3();
        const Body_applyLocalImpulse_relativePoint = new Vec3();
        const Body_updateMassProperties_halfExtents = new Vec3();

        /**
         * Base class for broadphase implementations
         * @author schteppe
         */
        class Broadphase {
          /**
           * The world to search for collisions in.
           */

          /**
           * If set to true, the broadphase uses bounding boxes for intersection tests, else it uses bounding spheres.
           */

          /**
           * Set to true if the objects in the world moved.
           */
          constructor() {
            this.world = void 0;
            this.useBoundingBoxes = void 0;
            this.dirty = void 0;
            this.world = null;
            this.useBoundingBoxes = false;
            this.dirty = true;
          }
          /**
           * Get the collision pairs from the world
           * @param world The world to search in
           * @param p1 Empty array to be filled with body objects
           * @param p2 Empty array to be filled with body objects
           */

          collisionPairs(world, p1, p2) {
            throw new Error(
              "collisionPairs not implemented for this BroadPhase class!"
            );
          }
          /**
           * Check if a body pair needs to be intersection tested at all.
           */

          needBroadphaseCollision(bodyA, bodyB) {
            // Check collision filter masks
            if (
              (bodyA.collisionFilterGroup & bodyB.collisionFilterMask) === 0 ||
              (bodyB.collisionFilterGroup & bodyA.collisionFilterMask) === 0
            ) {
              return false;
            } // Check types

            if (
              ((bodyA.type & Body.STATIC) !== 0 ||
                bodyA.sleepState === Body.SLEEPING) &&
              ((bodyB.type & Body.STATIC) !== 0 ||
                bodyB.sleepState === Body.SLEEPING)
            ) {
              // Both bodies are static or sleeping. Skip.
              return false;
            }

            return true;
          }
          /**
           * Check if the bounding volumes of two bodies intersect.
           */

          intersectionTest(bodyA, bodyB, pairs1, pairs2) {
            if (this.useBoundingBoxes) {
              this.doBoundingBoxBroadphase(bodyA, bodyB, pairs1, pairs2);
            } else {
              this.doBoundingSphereBroadphase(bodyA, bodyB, pairs1, pairs2);
            }
          }
          /**
           * Check if the bounding spheres of two bodies are intersecting.
           * @param pairs1 bodyA is appended to this array if intersection
           * @param pairs2 bodyB is appended to this array if intersection
           */

          doBoundingSphereBroadphase(bodyA, bodyB, pairs1, pairs2) {
            const r = Broadphase_collisionPairs_r;
            bodyB.position.vsub(bodyA.position, r);
            const boundingRadiusSum2 =
              (bodyA.boundingRadius + bodyB.boundingRadius) ** 2;
            const norm2 = r.lengthSquared();

            if (norm2 < boundingRadiusSum2) {
              pairs1.push(bodyA);
              pairs2.push(bodyB);
            }
          }
          /**
           * Check if the bounding boxes of two bodies are intersecting.
           */

          doBoundingBoxBroadphase(bodyA, bodyB, pairs1, pairs2) {
            if (bodyA.aabbNeedsUpdate) {
              bodyA.updateAABB();
            }

            if (bodyB.aabbNeedsUpdate) {
              bodyB.updateAABB();
            } // Check AABB / AABB

            if (bodyA.aabb.overlaps(bodyB.aabb)) {
              pairs1.push(bodyA);
              pairs2.push(bodyB);
            }
          }
          /**
           * Removes duplicate pairs from the pair arrays.
           */

          makePairsUnique(pairs1, pairs2) {
            const t = Broadphase_makePairsUnique_temp;
            const p1 = Broadphase_makePairsUnique_p1;
            const p2 = Broadphase_makePairsUnique_p2;
            const N = pairs1.length;

            for (let i = 0; i !== N; i++) {
              p1[i] = pairs1[i];
              p2[i] = pairs2[i];
            }

            pairs1.length = 0;
            pairs2.length = 0;

            for (let i = 0; i !== N; i++) {
              const id1 = p1[i].id;
              const id2 = p2[i].id;
              const key = id1 < id2 ? id1 + "," + id2 : id2 + "," + id1;
              t[key] = i;
              t.keys.push(key);
            }

            for (let i = 0; i !== t.keys.length; i++) {
              const key = t.keys.pop();
              const pairIndex = t[key];
              pairs1.push(p1[pairIndex]);
              pairs2.push(p2[pairIndex]);
              delete t[key];
            }
          }
          /**
           * To be implemented by subcasses
           */

          setWorld(world) {}
          /**
           * Check if the bounding spheres of two bodies overlap.
           */

          static boundingSphereCheck(bodyA, bodyB) {
            const dist = new Vec3(); // bsc_dist;

            bodyA.position.vsub(bodyB.position, dist);
            const sa = bodyA.shapes[0];
            const sb = bodyB.shapes[0];
            return (
              Math.pow(sa.boundingSphereRadius + sb.boundingSphereRadius, 2) >
              dist.lengthSquared()
            );
          }
          /**
           * Returns all the bodies within the AABB.
           */

          aabbQuery(world, aabb, result) {
            console.warn(
              ".aabbQuery is not implemented in this Broadphase subclass."
            );
            return [];
          }
        } // Temp objects

        const Broadphase_collisionPairs_r = new Vec3();
        const Broadphase_makePairsUnique_temp = {
          keys: [],
        };
        const Broadphase_makePairsUnique_p1 = [];
        const Broadphase_makePairsUnique_p2 = [];

        /**
         * Axis aligned uniform grid broadphase.
         * @todo Needs support for more than just planes and spheres.
         */
        class GridBroadphase extends Broadphase {
          /**
           * Number of boxes along x
           */

          /**
           * Number of boxes along y
           */

          /**
           * Number of boxes along z
           */

          /**
           * aabbMin
           */

          /**
           * aabbMax
           */

          /**
           * bins
           */

          /**
           * binLengths
           */

          /**
           * @param nx Number of boxes along x.
           * @param ny Number of boxes along y.
           * @param nz Number of boxes along z.
           */
          constructor(
            aabbMin = new Vec3(100, 100, 100),
            aabbMax = new Vec3(-100, -100, -100),
            nx = 10,
            ny = 10,
            nz = 10
          ) {
            super();
            this.nx = void 0;
            this.ny = void 0;
            this.nz = void 0;
            this.aabbMin = void 0;
            this.aabbMax = void 0;
            this.bins = void 0;
            this.binLengths = void 0;
            this.nx = nx;
            this.ny = ny;
            this.nz = nz;
            this.aabbMin = aabbMin;
            this.aabbMax = aabbMax;
            const nbins = this.nx * this.ny * this.nz;

            if (nbins <= 0) {
              throw "GridBroadphase: Each dimension's n must be >0";
            }

            this.bins = [];
            this.binLengths = []; // Rather than continually resizing arrays (thrashing the memory), just record length and allow them to grow

            this.bins.length = nbins;
            this.binLengths.length = nbins;

            for (let i = 0; i < nbins; i++) {
              this.bins[i] = [];
              this.binLengths[i] = 0;
            }
          }
          /**
           * Get all the collision pairs in the physics world
           */

          collisionPairs(world, pairs1, pairs2) {
            const N = world.numObjects();
            const bodies = world.bodies;
            const max = this.aabbMax;
            const min = this.aabbMin;
            const nx = this.nx;
            const ny = this.ny;
            const nz = this.nz;
            const xstep = ny * nz;
            const ystep = nz;
            const zstep = 1;
            const xmax = max.x;
            const ymax = max.y;
            const zmax = max.z;
            const xmin = min.x;
            const ymin = min.y;
            const zmin = min.z;
            const xmult = nx / (xmax - xmin);
            const ymult = ny / (ymax - ymin);
            const zmult = nz / (zmax - zmin);
            const binsizeX = (xmax - xmin) / nx;
            const binsizeY = (ymax - ymin) / ny;
            const binsizeZ = (zmax - zmin) / nz;
            const binRadius =
              Math.sqrt(
                binsizeX * binsizeX + binsizeY * binsizeY + binsizeZ * binsizeZ
              ) * 0.5;
            const types = Shape.types;
            const SPHERE = types.SPHERE;
            const PLANE = types.PLANE;
            types.BOX;
            types.COMPOUND;
            types.CONVEXPOLYHEDRON;
            const bins = this.bins;
            const binLengths = this.binLengths;
            const Nbins = this.bins.length; // Reset bins

            for (let i = 0; i !== Nbins; i++) {
              binLengths[i] = 0;
            }

            const ceil = Math.ceil;

            function addBoxToBins(x0, y0, z0, x1, y1, z1, bi) {
              let xoff0 = ((x0 - xmin) * xmult) | 0;
              let yoff0 = ((y0 - ymin) * ymult) | 0;
              let zoff0 = ((z0 - zmin) * zmult) | 0;
              let xoff1 = ceil((x1 - xmin) * xmult);
              let yoff1 = ceil((y1 - ymin) * ymult);
              let zoff1 = ceil((z1 - zmin) * zmult);

              if (xoff0 < 0) {
                xoff0 = 0;
              } else if (xoff0 >= nx) {
                xoff0 = nx - 1;
              }

              if (yoff0 < 0) {
                yoff0 = 0;
              } else if (yoff0 >= ny) {
                yoff0 = ny - 1;
              }

              if (zoff0 < 0) {
                zoff0 = 0;
              } else if (zoff0 >= nz) {
                zoff0 = nz - 1;
              }

              if (xoff1 < 0) {
                xoff1 = 0;
              } else if (xoff1 >= nx) {
                xoff1 = nx - 1;
              }

              if (yoff1 < 0) {
                yoff1 = 0;
              } else if (yoff1 >= ny) {
                yoff1 = ny - 1;
              }

              if (zoff1 < 0) {
                zoff1 = 0;
              } else if (zoff1 >= nz) {
                zoff1 = nz - 1;
              }

              xoff0 *= xstep;
              yoff0 *= ystep;
              zoff0 *= zstep;
              xoff1 *= xstep;
              yoff1 *= ystep;
              zoff1 *= zstep;

              for (let xoff = xoff0; xoff <= xoff1; xoff += xstep) {
                for (let yoff = yoff0; yoff <= yoff1; yoff += ystep) {
                  for (let zoff = zoff0; zoff <= zoff1; zoff += zstep) {
                    const idx = xoff + yoff + zoff;
                    bins[idx][binLengths[idx]++] = bi;
                  }
                }
              }
            } // Put all bodies into the bins

            for (let i = 0; i !== N; i++) {
              const bi = bodies[i];
              const si = bi.shapes[0];

              switch (si.type) {
                case SPHERE: {
                  const shape = si; // Put in bin
                  // check if overlap with other bins

                  const x = bi.position.x;
                  const y = bi.position.y;
                  const z = bi.position.z;
                  const r = shape.radius;
                  addBoxToBins(x - r, y - r, z - r, x + r, y + r, z + r, bi);
                  break;
                }

                case PLANE: {
                  const shape = si;

                  if (shape.worldNormalNeedsUpdate) {
                    shape.computeWorldNormal(bi.quaternion);
                  }

                  const planeNormal = shape.worldNormal; //Relative position from origin of plane object to the first bin
                  //Incremented as we iterate through the bins

                  const xreset = xmin + binsizeX * 0.5 - bi.position.x;
                  const yreset = ymin + binsizeY * 0.5 - bi.position.y;
                  const zreset = zmin + binsizeZ * 0.5 - bi.position.z;
                  const d = GridBroadphase_collisionPairs_d;
                  d.set(xreset, yreset, zreset);

                  for (
                    let xi = 0, xoff = 0;
                    xi !== nx;
                    xi++, xoff += xstep, d.y = yreset, d.x += binsizeX
                  ) {
                    for (
                      let yi = 0, yoff = 0;
                      yi !== ny;
                      yi++, yoff += ystep, d.z = zreset, d.y += binsizeY
                    ) {
                      for (
                        let zi = 0, zoff = 0;
                        zi !== nz;
                        zi++, zoff += zstep, d.z += binsizeZ
                      ) {
                        if (d.dot(planeNormal) < binRadius) {
                          const idx = xoff + yoff + zoff;
                          bins[idx][binLengths[idx]++] = bi;
                        }
                      }
                    }
                  }

                  break;
                }

                default: {
                  if (bi.aabbNeedsUpdate) {
                    bi.updateAABB();
                  }

                  addBoxToBins(
                    bi.aabb.lowerBound.x,
                    bi.aabb.lowerBound.y,
                    bi.aabb.lowerBound.z,
                    bi.aabb.upperBound.x,
                    bi.aabb.upperBound.y,
                    bi.aabb.upperBound.z,
                    bi
                  );
                  break;
                }
              }
            } // Check each bin

            for (let i = 0; i !== Nbins; i++) {
              const binLength = binLengths[i]; //Skip bins with no potential collisions

              if (binLength > 1) {
                const bin = bins[i]; // Do N^2 broadphase inside

                for (let xi = 0; xi !== binLength; xi++) {
                  const bi = bin[xi];

                  for (let yi = 0; yi !== xi; yi++) {
                    const bj = bin[yi];

                    if (this.needBroadphaseCollision(bi, bj)) {
                      this.intersectionTest(bi, bj, pairs1, pairs2);
                    }
                  }
                }
              }
            } //	for (let zi = 0, zoff=0; zi < nz; zi++, zoff+= zstep) {
            //		console.log("layer "+zi);
            //		for (let yi = 0, yoff=0; yi < ny; yi++, yoff += ystep) {
            //			const row = '';
            //			for (let xi = 0, xoff=0; xi < nx; xi++, xoff += xstep) {
            //				const idx = xoff + yoff + zoff;
            //				row += ' ' + binLengths[idx];
            //			}
            //			console.log(row);
            //		}
            //	}

            this.makePairsUnique(pairs1, pairs2);
          }
        }
        const GridBroadphase_collisionPairs_d = new Vec3();

        /**
         * Naive broadphase implementation, used in lack of better ones.
         *
         * The naive broadphase looks at all possible pairs without restriction, therefore it has complexity N^2 _(which is bad)_
         */
        class NaiveBroadphase extends Broadphase {
          /**
           * @todo Remove useless constructor
           */
          constructor() {
            super();
          }
          /**
           * Get all the collision pairs in the physics world
           */

          collisionPairs(world, pairs1, pairs2) {
            const bodies = world.bodies;
            const n = bodies.length;
            let bi;
            let bj; // Naive N^2 ftw!

            for (let i = 0; i !== n; i++) {
              for (let j = 0; j !== i; j++) {
                bi = bodies[i];
                bj = bodies[j];

                if (!this.needBroadphaseCollision(bi, bj)) {
                  continue;
                }

                this.intersectionTest(bi, bj, pairs1, pairs2);
              }
            }
          }
          /**
           * Returns all the bodies within an AABB.
           * @param result An array to store resulting bodies in.
           */

          aabbQuery(world, aabb, result = []) {
            for (let i = 0; i < world.bodies.length; i++) {
              const b = world.bodies[i];

              if (b.aabbNeedsUpdate) {
                b.updateAABB();
              } // Ugly hack until Body gets aabb

              if (b.aabb.overlaps(aabb)) {
                result.push(b);
              }
            }

            return result;
          }
        }

        /**
         * Storage for Ray casting data
         */
        class RaycastResult {
          /**
           * rayFromWorld
           */

          /**
           * rayToWorld
           */

          /**
           * hitNormalWorld
           */

          /**
           * hitPointWorld
           */

          /**
           * hasHit
           */

          /**
           * shape
           */

          /**
           * body
           */

          /**
           * The index of the hit triangle, if the hit shape was a trimesh
           */

          /**
           * Distance to the hit. Will be set to -1 if there was no hit
           */

          /**
           * If the ray should stop traversing the bodies
           */
          constructor() {
            this.rayFromWorld = void 0;
            this.rayToWorld = void 0;
            this.hitNormalWorld = void 0;
            this.hitPointWorld = void 0;
            this.hasHit = void 0;
            this.shape = void 0;
            this.body = void 0;
            this.hitFaceIndex = void 0;
            this.distance = void 0;
            this.shouldStop = void 0;
            this.rayFromWorld = new Vec3();
            this.rayToWorld = new Vec3();
            this.hitNormalWorld = new Vec3();
            this.hitPointWorld = new Vec3();
            this.hasHit = false;
            this.shape = null;
            this.body = null;
            this.hitFaceIndex = -1;
            this.distance = -1;
            this.shouldStop = false;
          }
          /**
           * Reset all result data.
           */

          reset() {
            this.rayFromWorld.setZero();
            this.rayToWorld.setZero();
            this.hitNormalWorld.setZero();
            this.hitPointWorld.setZero();
            this.hasHit = false;
            this.shape = null;
            this.body = null;
            this.hitFaceIndex = -1;
            this.distance = -1;
            this.shouldStop = false;
          }
          /**
           * abort
           */

          abort() {
            this.shouldStop = true;
          }
          /**
           * Set result data.
           */

          set(
            rayFromWorld,
            rayToWorld,
            hitNormalWorld,
            hitPointWorld,
            shape,
            body,
            distance
          ) {
            this.rayFromWorld.copy(rayFromWorld);
            this.rayToWorld.copy(rayToWorld);
            this.hitNormalWorld.copy(hitNormalWorld);
            this.hitPointWorld.copy(hitPointWorld);
            this.shape = shape;
            this.body = body;
            this.distance = distance;
          }
        }

        let _Shape$types$SPHERE,
          _Shape$types$PLANE,
          _Shape$types$BOX,
          _Shape$types$CYLINDER,
          _Shape$types$CONVEXPO,
          _Shape$types$HEIGHTFI,
          _Shape$types$TRIMESH;

        /**
         * RAY_MODES
         */
        const RAY_MODES = {
          /** CLOSEST */
          CLOSEST: 1,

          /** ANY */
          ANY: 2,

          /** ALL */
          ALL: 4,
        };
        /**
         * RayMode
         */

        _Shape$types$SPHERE = Shape.types.SPHERE;
        _Shape$types$PLANE = Shape.types.PLANE;
        _Shape$types$BOX = Shape.types.BOX;
        _Shape$types$CYLINDER = Shape.types.CYLINDER;
        _Shape$types$CONVEXPO = Shape.types.CONVEXPOLYHEDRON;
        _Shape$types$HEIGHTFI = Shape.types.HEIGHTFIELD;
        _Shape$types$TRIMESH = Shape.types.TRIMESH;

        /**
         * A line in 3D space that intersects bodies and return points.
         */
        class Ray {
          /**
           * from
           */

          /**
           * to
           */

          /**
           * direction
           */

          /**
           * The precision of the ray. Used when checking parallelity etc.
           * @default 0.0001
           */

          /**
           * Set to `false` if you don't want the Ray to take `collisionResponse` flags into account on bodies and shapes.
           * @default true
           */

          /**
           * If set to `true`, the ray skips any hits with normal.dot(rayDirection) < 0.
           * @default false
           */

          /**
           * collisionFilterMask
           * @default -1
           */

          /**
           * collisionFilterGroup
           * @default -1
           */

          /**
           * The intersection mode. Should be Ray.ANY, Ray.ALL or Ray.CLOSEST.
           * @default RAY.ANY
           */

          /**
           * Current result object.
           */

          /**
           * Will be set to `true` during intersectWorld() if the ray hit anything.
           */

          /**
           * User-provided result callback. Will be used if mode is Ray.ALL.
           */

          /**
           * CLOSEST
           */

          /**
           * ANY
           */

          /**
           * ALL
           */
          get [_Shape$types$SPHERE]() {
            return this._intersectSphere;
          }

          get [_Shape$types$PLANE]() {
            return this._intersectPlane;
          }

          get [_Shape$types$BOX]() {
            return this._intersectBox;
          }

          get [_Shape$types$CYLINDER]() {
            return this._intersectConvex;
          }

          get [_Shape$types$CONVEXPO]() {
            return this._intersectConvex;
          }

          get [_Shape$types$HEIGHTFI]() {
            return this._intersectHeightfield;
          }

          get [_Shape$types$TRIMESH]() {
            return this._intersectTrimesh;
          }

          constructor(from = new Vec3(), to = new Vec3()) {
            this.from = void 0;
            this.to = void 0;
            this.direction = void 0;
            this.precision = void 0;
            this.checkCollisionResponse = void 0;
            this.skipBackfaces = void 0;
            this.collisionFilterMask = void 0;
            this.collisionFilterGroup = void 0;
            this.mode = void 0;
            this.result = void 0;
            this.hasHit = void 0;
            this.callback = void 0;
            this.from = from.clone();
            this.to = to.clone();
            this.direction = new Vec3();
            this.precision = 0.0001;
            this.checkCollisionResponse = true;
            this.skipBackfaces = false;
            this.collisionFilterMask = -1;
            this.collisionFilterGroup = -1;
            this.mode = Ray.ANY;
            this.result = new RaycastResult();
            this.hasHit = false;

            this.callback = (result) => {};
          }
          /**
           * Do itersection against all bodies in the given World.
           * @return True if the ray hit anything, otherwise false.
           */

          intersectWorld(world, options) {
            this.mode = options.mode || Ray.ANY;
            this.result = options.result || new RaycastResult();
            this.skipBackfaces = !!options.skipBackfaces;
            this.collisionFilterMask =
              typeof options.collisionFilterMask !== "undefined"
                ? options.collisionFilterMask
                : -1;
            this.collisionFilterGroup =
              typeof options.collisionFilterGroup !== "undefined"
                ? options.collisionFilterGroup
                : -1;
            this.checkCollisionResponse =
              typeof options.checkCollisionResponse !== "undefined"
                ? options.checkCollisionResponse
                : true;

            if (options.from) {
              this.from.copy(options.from);
            }

            if (options.to) {
              this.to.copy(options.to);
            }

            this.callback = options.callback || (() => {});

            this.hasHit = false;
            this.result.reset();
            this.updateDirection();
            this.getAABB(tmpAABB$1);
            tmpArray.length = 0;
            world.broadphase.aabbQuery(world, tmpAABB$1, tmpArray);
            this.intersectBodies(tmpArray);
            return this.hasHit;
          }
          /**
           * Shoot a ray at a body, get back information about the hit.
           * @deprecated @param result set the result property of the Ray instead.
           */

          intersectBody(body, result) {
            if (result) {
              this.result = result;
              this.updateDirection();
            }

            const checkCollisionResponse = this.checkCollisionResponse;

            if (checkCollisionResponse && !body.collisionResponse) {
              return;
            }

            if (
              (this.collisionFilterGroup & body.collisionFilterMask) === 0 ||
              (body.collisionFilterGroup & this.collisionFilterMask) === 0
            ) {
              return;
            }

            const xi = intersectBody_xi;
            const qi = intersectBody_qi;

            for (let i = 0, N = body.shapes.length; i < N; i++) {
              const shape = body.shapes[i];

              if (checkCollisionResponse && !shape.collisionResponse) {
                continue; // Skip
              }

              body.quaternion.mult(body.shapeOrientations[i], qi);
              body.quaternion.vmult(body.shapeOffsets[i], xi);
              xi.vadd(body.position, xi);
              this.intersectShape(shape, qi, xi, body);

              if (this.result.shouldStop) {
                break;
              }
            }
          }
          /**
           * Shoot a ray at an array bodies, get back information about the hit.
           * @param bodies An array of Body objects.
           * @deprecated @param result set the result property of the Ray instead.
           *
           */

          intersectBodies(bodies, result) {
            if (result) {
              this.result = result;
              this.updateDirection();
            }

            for (
              let i = 0, l = bodies.length;
              !this.result.shouldStop && i < l;
              i++
            ) {
              this.intersectBody(bodies[i]);
            }
          }
          /**
           * Updates the direction vector.
           */

          updateDirection() {
            this.to.vsub(this.from, this.direction);
            this.direction.normalize();
          }

          intersectShape(shape, quat, position, body) {
            const from = this.from; // Checking boundingSphere

            const distance = distanceFromIntersection(
              from,
              this.direction,
              position
            );

            if (distance > shape.boundingSphereRadius) {
              return;
            }

            const intersectMethod = this[shape.type];

            if (intersectMethod) {
              intersectMethod.call(this, shape, quat, position, body, shape);
            }
          }

          _intersectBox(box, quat, position, body, reportedShape) {
            return this._intersectConvex(
              box.convexPolyhedronRepresentation,
              quat,
              position,
              body,
              reportedShape
            );
          }

          _intersectPlane(shape, quat, position, body, reportedShape) {
            const from = this.from;
            const to = this.to;
            const direction = this.direction; // Get plane normal

            const worldNormal = new Vec3(0, 0, 1);
            quat.vmult(worldNormal, worldNormal);
            const len = new Vec3();
            from.vsub(position, len);
            const planeToFrom = len.dot(worldNormal);
            to.vsub(position, len);
            const planeToTo = len.dot(worldNormal);

            if (planeToFrom * planeToTo > 0) {
              // "from" and "to" are on the same side of the plane... bail out
              return;
            }

            if (from.distanceTo(to) < planeToFrom) {
              return;
            }

            const n_dot_dir = worldNormal.dot(direction);

            if (Math.abs(n_dot_dir) < this.precision) {
              // No intersection
              return;
            }

            const planePointToFrom = new Vec3();
            const dir_scaled_with_t = new Vec3();
            const hitPointWorld = new Vec3();
            from.vsub(position, planePointToFrom);
            const t = -worldNormal.dot(planePointToFrom) / n_dot_dir;
            direction.scale(t, dir_scaled_with_t);
            from.vadd(dir_scaled_with_t, hitPointWorld);
            this.reportIntersection(
              worldNormal,
              hitPointWorld,
              reportedShape,
              body,
              -1
            );
          }
          /**
           * Get the world AABB of the ray.
           */

          getAABB(aabb) {
            const { lowerBound, upperBound } = aabb;
            const to = this.to;
            const from = this.from;
            lowerBound.x = Math.min(to.x, from.x);
            lowerBound.y = Math.min(to.y, from.y);
            lowerBound.z = Math.min(to.z, from.z);
            upperBound.x = Math.max(to.x, from.x);
            upperBound.y = Math.max(to.y, from.y);
            upperBound.z = Math.max(to.z, from.z);
          }

          _intersectHeightfield(shape, quat, position, body, reportedShape) {
            shape.data;
            shape.elementSize; // Convert the ray to local heightfield coordinates

            const localRay = intersectHeightfield_localRay; //new Ray(this.from, this.to);

            localRay.from.copy(this.from);
            localRay.to.copy(this.to);
            Transform.pointToLocalFrame(
              position,
              quat,
              localRay.from,
              localRay.from
            );
            Transform.pointToLocalFrame(
              position,
              quat,
              localRay.to,
              localRay.to
            );
            localRay.updateDirection(); // Get the index of the data points to test against

            const index = intersectHeightfield_index;
            let iMinX;
            let iMinY;
            let iMaxX;
            let iMaxY; // Set to max

            iMinX = iMinY = 0;
            iMaxX = iMaxY = shape.data.length - 1;
            const aabb = new AABB();
            localRay.getAABB(aabb);
            shape.getIndexOfPosition(
              aabb.lowerBound.x,
              aabb.lowerBound.y,
              index,
              true
            );
            iMinX = Math.max(iMinX, index[0]);
            iMinY = Math.max(iMinY, index[1]);
            shape.getIndexOfPosition(
              aabb.upperBound.x,
              aabb.upperBound.y,
              index,
              true
            );
            iMaxX = Math.min(iMaxX, index[0] + 1);
            iMaxY = Math.min(iMaxY, index[1] + 1);

            for (let i = iMinX; i < iMaxX; i++) {
              for (let j = iMinY; j < iMaxY; j++) {
                if (this.result.shouldStop) {
                  return;
                }

                shape.getAabbAtIndex(i, j, aabb);

                if (!aabb.overlapsRay(localRay)) {
                  continue;
                } // Lower triangle

                shape.getConvexTrianglePillar(i, j, false);
                Transform.pointToWorldFrame(
                  position,
                  quat,
                  shape.pillarOffset,
                  worldPillarOffset
                );

                this._intersectConvex(
                  shape.pillarConvex,
                  quat,
                  worldPillarOffset,
                  body,
                  reportedShape,
                  intersectConvexOptions
                );

                if (this.result.shouldStop) {
                  return;
                } // Upper triangle

                shape.getConvexTrianglePillar(i, j, true);
                Transform.pointToWorldFrame(
                  position,
                  quat,
                  shape.pillarOffset,
                  worldPillarOffset
                );

                this._intersectConvex(
                  shape.pillarConvex,
                  quat,
                  worldPillarOffset,
                  body,
                  reportedShape,
                  intersectConvexOptions
                );
              }
            }
          }

          _intersectSphere(sphere, quat, position, body, reportedShape) {
            const from = this.from;
            const to = this.to;
            const r = sphere.radius;
            const a =
              (to.x - from.x) ** 2 +
              (to.y - from.y) ** 2 +
              (to.z - from.z) ** 2;
            const b =
              2 *
              ((to.x - from.x) * (from.x - position.x) +
                (to.y - from.y) * (from.y - position.y) +
                (to.z - from.z) * (from.z - position.z));
            const c =
              (from.x - position.x) ** 2 +
              (from.y - position.y) ** 2 +
              (from.z - position.z) ** 2 -
              r ** 2;
            const delta = b ** 2 - 4 * a * c;
            const intersectionPoint = Ray_intersectSphere_intersectionPoint;
            const normal = Ray_intersectSphere_normal;

            if (delta < 0) {
              // No intersection
              return;
            } else if (delta === 0) {
              // single intersection point
              from.lerp(to, delta, intersectionPoint);
              intersectionPoint.vsub(position, normal);
              normal.normalize();
              this.reportIntersection(
                normal,
                intersectionPoint,
                reportedShape,
                body,
                -1
              );
            } else {
              const d1 = (-b - Math.sqrt(delta)) / (2 * a);
              const d2 = (-b + Math.sqrt(delta)) / (2 * a);

              if (d1 >= 0 && d1 <= 1) {
                from.lerp(to, d1, intersectionPoint);
                intersectionPoint.vsub(position, normal);
                normal.normalize();
                this.reportIntersection(
                  normal,
                  intersectionPoint,
                  reportedShape,
                  body,
                  -1
                );
              }

              if (this.result.shouldStop) {
                return;
              }

              if (d2 >= 0 && d2 <= 1) {
                from.lerp(to, d2, intersectionPoint);
                intersectionPoint.vsub(position, normal);
                normal.normalize();
                this.reportIntersection(
                  normal,
                  intersectionPoint,
                  reportedShape,
                  body,
                  -1
                );
              }
            }
          }

          _intersectConvex(
            shape,
            quat,
            position,
            body,
            reportedShape,
            options
          ) {
            const normal = intersectConvex_normal;
            const vector = intersectConvex_vector;
            const faceList = (options && options.faceList) || null; // Checking faces

            const faces = shape.faces;
            const vertices = shape.vertices;
            const normals = shape.faceNormals;
            const direction = this.direction;
            const from = this.from;
            const to = this.to;
            const fromToDistance = from.distanceTo(to);
            const Nfaces = faceList ? faceList.length : faces.length;
            const result = this.result;

            for (let j = 0; !result.shouldStop && j < Nfaces; j++) {
              const fi = faceList ? faceList[j] : j;
              const face = faces[fi];
              const faceNormal = normals[fi];
              const q = quat;
              const x = position; // determine if ray intersects the plane of the face
              // note: this works regardless of the direction of the face normal
              // Get plane point in world coordinates...

              vector.copy(vertices[face[0]]);
              q.vmult(vector, vector);
              vector.vadd(x, vector); // ...but make it relative to the ray from. We'll fix this later.

              vector.vsub(from, vector); // Get plane normal

              q.vmult(faceNormal, normal); // If this dot product is negative, we have something interesting

              const dot = direction.dot(normal); // Bail out if ray and plane are parallel

              if (Math.abs(dot) < this.precision) {
                continue;
              } // calc distance to plane

              const scalar = normal.dot(vector) / dot; // if negative distance, then plane is behind ray

              if (scalar < 0) {
                continue;
              } // if (dot < 0) {
              // Intersection point is from + direction * scalar

              direction.scale(scalar, intersectPoint);
              intersectPoint.vadd(from, intersectPoint); // a is the point we compare points b and c with.

              a.copy(vertices[face[0]]);
              q.vmult(a, a);
              x.vadd(a, a);

              for (let i = 1; !result.shouldStop && i < face.length - 1; i++) {
                // Transform 3 vertices to world coords
                b.copy(vertices[face[i]]);
                c.copy(vertices[face[i + 1]]);
                q.vmult(b, b);
                q.vmult(c, c);
                x.vadd(b, b);
                x.vadd(c, c);
                const distance = intersectPoint.distanceTo(from);

                if (
                  !(
                    Ray.pointInTriangle(intersectPoint, a, b, c) ||
                    Ray.pointInTriangle(intersectPoint, b, a, c)
                  ) ||
                  distance > fromToDistance
                ) {
                  continue;
                }

                this.reportIntersection(
                  normal,
                  intersectPoint,
                  reportedShape,
                  body,
                  fi
                );
              } // }
            }
          }
          /**
           * @todo Optimize by transforming the world to local space first.
           * @todo Use Octree lookup
           */

          _intersectTrimesh(
            mesh,
            quat,
            position,
            body,
            reportedShape,
            options
          ) {
            const normal = intersectTrimesh_normal;
            const triangles = intersectTrimesh_triangles;
            const treeTransform = intersectTrimesh_treeTransform;
            const vector = intersectConvex_vector;
            const localDirection = intersectTrimesh_localDirection;
            const localFrom = intersectTrimesh_localFrom;
            const localTo = intersectTrimesh_localTo;
            const worldIntersectPoint = intersectTrimesh_worldIntersectPoint;
            const worldNormal = intersectTrimesh_worldNormal; // Checking faces

            const indices = mesh.indices;
            mesh.vertices; // const normals = mesh.faceNormals

            const from = this.from;
            const to = this.to;
            const direction = this.direction;
            treeTransform.position.copy(position);
            treeTransform.quaternion.copy(quat); // Transform ray to local space!

            Transform.vectorToLocalFrame(
              position,
              quat,
              direction,
              localDirection
            );
            Transform.pointToLocalFrame(position, quat, from, localFrom);
            Transform.pointToLocalFrame(position, quat, to, localTo);
            localTo.x *= mesh.scale.x;
            localTo.y *= mesh.scale.y;
            localTo.z *= mesh.scale.z;
            localFrom.x *= mesh.scale.x;
            localFrom.y *= mesh.scale.y;
            localFrom.z *= mesh.scale.z;
            localTo.vsub(localFrom, localDirection);
            localDirection.normalize();
            const fromToDistanceSquared = localFrom.distanceSquared(localTo);
            mesh.tree.rayQuery(this, treeTransform, triangles);

            for (
              let i = 0, N = triangles.length;
              !this.result.shouldStop && i !== N;
              i++
            ) {
              const trianglesIndex = triangles[i];
              mesh.getNormal(trianglesIndex, normal); // determine if ray intersects the plane of the face
              // note: this works regardless of the direction of the face normal
              // Get plane point in world coordinates...

              mesh.getVertex(indices[trianglesIndex * 3], a); // ...but make it relative to the ray from. We'll fix this later.

              a.vsub(localFrom, vector); // If this dot product is negative, we have something interesting

              const dot = localDirection.dot(normal); // Bail out if ray and plane are parallel
              // if (Math.abs( dot ) < this.precision){
              //     continue;
              // }
              // calc distance to plane

              const scalar = normal.dot(vector) / dot; // if negative distance, then plane is behind ray

              if (scalar < 0) {
                continue;
              } // Intersection point is from + direction * scalar

              localDirection.scale(scalar, intersectPoint);
              intersectPoint.vadd(localFrom, intersectPoint); // Get triangle vertices

              mesh.getVertex(indices[trianglesIndex * 3 + 1], b);
              mesh.getVertex(indices[trianglesIndex * 3 + 2], c);
              const squaredDistance = intersectPoint.distanceSquared(localFrom);

              if (
                !(
                  Ray.pointInTriangle(intersectPoint, b, a, c) ||
                  Ray.pointInTriangle(intersectPoint, a, b, c)
                ) ||
                squaredDistance > fromToDistanceSquared
              ) {
                continue;
              } // transform intersectpoint and normal to world

              Transform.vectorToWorldFrame(quat, normal, worldNormal);
              Transform.pointToWorldFrame(
                position,
                quat,
                intersectPoint,
                worldIntersectPoint
              );
              this.reportIntersection(
                worldNormal,
                worldIntersectPoint,
                reportedShape,
                body,
                trianglesIndex
              );
            }

            triangles.length = 0;
          }
          /**
           * @return True if the intersections should continue
           */

          reportIntersection(normal, hitPointWorld, shape, body, hitFaceIndex) {
            const from = this.from;
            const to = this.to;
            const distance = from.distanceTo(hitPointWorld);
            const result = this.result; // Skip back faces?

            if (this.skipBackfaces && normal.dot(this.direction) > 0) {
              return;
            }

            result.hitFaceIndex =
              typeof hitFaceIndex !== "undefined" ? hitFaceIndex : -1;

            switch (this.mode) {
              case Ray.ALL:
                this.hasHit = true;
                result.set(
                  from,
                  to,
                  normal,
                  hitPointWorld,
                  shape,
                  body,
                  distance
                );
                result.hasHit = true;
                this.callback(result);
                break;

              case Ray.CLOSEST:
                // Store if closer than current closest
                if (distance < result.distance || !result.hasHit) {
                  this.hasHit = true;
                  result.hasHit = true;
                  result.set(
                    from,
                    to,
                    normal,
                    hitPointWorld,
                    shape,
                    body,
                    distance
                  );
                }

                break;

              case Ray.ANY:
                // Report and stop.
                this.hasHit = true;
                result.hasHit = true;
                result.set(
                  from,
                  to,
                  normal,
                  hitPointWorld,
                  shape,
                  body,
                  distance
                );
                result.shouldStop = true;
                break;
            }
          }
          /**
           * As per "Barycentric Technique" as named
           * {@link https://www.blackpawn.com/texts/pointinpoly/default.html here} but without the division
           */

          static pointInTriangle(p, a, b, c) {
            c.vsub(a, v0);
            b.vsub(a, v1);
            p.vsub(a, v2);
            const dot00 = v0.dot(v0);
            const dot01 = v0.dot(v1);
            const dot02 = v0.dot(v2);
            const dot11 = v1.dot(v1);
            const dot12 = v1.dot(v2);
            let u;
            let v;
            return (
              (u = dot11 * dot02 - dot01 * dot12) >= 0 &&
              (v = dot00 * dot12 - dot01 * dot02) >= 0 &&
              u + v < dot00 * dot11 - dot01 * dot01
            );
          }
        }
        Ray.CLOSEST = RAY_MODES.CLOSEST;
        Ray.ANY = RAY_MODES.ANY;
        Ray.ALL = RAY_MODES.ALL;
        const tmpAABB$1 = new AABB();
        const tmpArray = [];
        const v1 = new Vec3();
        const v2 = new Vec3();
        const intersectBody_xi = new Vec3();
        const intersectBody_qi = new Quaternion();
        const intersectPoint = new Vec3();
        const a = new Vec3();
        const b = new Vec3();
        const c = new Vec3();
        const intersectConvexOptions = {
          faceList: [0],
        };
        const worldPillarOffset = new Vec3();
        const intersectHeightfield_localRay = new Ray();
        const intersectHeightfield_index = [];
        const Ray_intersectSphere_intersectionPoint = new Vec3();
        const Ray_intersectSphere_normal = new Vec3();
        const intersectConvex_normal = new Vec3();
        const intersectConvex_vector = new Vec3();
        const intersectTrimesh_normal = new Vec3();
        const intersectTrimesh_localDirection = new Vec3();
        const intersectTrimesh_localFrom = new Vec3();
        const intersectTrimesh_localTo = new Vec3();
        const intersectTrimesh_worldNormal = new Vec3();
        const intersectTrimesh_worldIntersectPoint = new Vec3();
        new AABB();
        const intersectTrimesh_triangles = [];
        const intersectTrimesh_treeTransform = new Transform();
        const v0 = new Vec3();
        const intersect = new Vec3();

        function distanceFromIntersection(from, direction, position) {
          // v0 is vector from from to position
          position.vsub(from, v0);
          const dot = v0.dot(direction); // intersect = direction*dot + from

          direction.scale(dot, intersect);
          intersect.vadd(from, intersect);
          const distance = position.distanceTo(intersect);
          return distance;
        }

        /**
         * Sweep and prune broadphase along one axis.
         */
        class SAPBroadphase extends Broadphase {
          /**
           * List of bodies currently in the broadphase.
           */

          /**
           * The world to search in.
           */

          /**
           * Axis to sort the bodies along.
           * Set to 0 for x axis, and 1 for y axis.
           * For best performance, pick the axis where bodies are most distributed.
           */

          /**
           * Check if the bounds of two bodies overlap, along the given SAP axis.
           */
          static checkBounds(bi, bj, axisIndex) {
            let biPos;
            let bjPos;

            if (axisIndex === 0) {
              biPos = bi.position.x;
              bjPos = bj.position.x;
            } else if (axisIndex === 1) {
              biPos = bi.position.y;
              bjPos = bj.position.y;
            } else if (axisIndex === 2) {
              biPos = bi.position.z;
              bjPos = bj.position.z;
            }

            const ri = bi.boundingRadius,
              rj = bj.boundingRadius,
              boundA2 = biPos + ri,
              boundB1 = bjPos - rj;
            return boundB1 < boundA2;
          } // Note: these are identical, save for x/y/z lowerbound

          /**
           * insertionSortX
           */

          static insertionSortX(a) {
            for (let i = 1, l = a.length; i < l; i++) {
              const v = a[i];
              let j;

              for (j = i - 1; j >= 0; j--) {
                if (a[j].aabb.lowerBound.x <= v.aabb.lowerBound.x) {
                  break;
                }

                a[j + 1] = a[j];
              }

              a[j + 1] = v;
            }

            return a;
          }
          /**
           * insertionSortY
           */

          static insertionSortY(a) {
            for (let i = 1, l = a.length; i < l; i++) {
              const v = a[i];
              let j;

              for (j = i - 1; j >= 0; j--) {
                if (a[j].aabb.lowerBound.y <= v.aabb.lowerBound.y) {
                  break;
                }

                a[j + 1] = a[j];
              }

              a[j + 1] = v;
            }

            return a;
          }
          /**
           * insertionSortZ
           */

          static insertionSortZ(a) {
            for (let i = 1, l = a.length; i < l; i++) {
              const v = a[i];
              let j;

              for (j = i - 1; j >= 0; j--) {
                if (a[j].aabb.lowerBound.z <= v.aabb.lowerBound.z) {
                  break;
                }

                a[j + 1] = a[j];
              }

              a[j + 1] = v;
            }

            return a;
          }

          constructor(world) {
            super();
            this.axisList = void 0;
            this.world = void 0;
            this.axisIndex = void 0;
            this._addBodyHandler = void 0;
            this._removeBodyHandler = void 0;
            this.axisList = [];
            this.world = null;
            this.axisIndex = 0;
            const axisList = this.axisList;

            this._addBodyHandler = (event) => {
              axisList.push(event.body);
            };

            this._removeBodyHandler = (event) => {
              const idx = axisList.indexOf(event.body);

              if (idx !== -1) {
                axisList.splice(idx, 1);
              }
            };

            if (world) {
              this.setWorld(world);
            }
          }
          /**
           * Change the world
           */

          setWorld(world) {
            // Clear the old axis array
            this.axisList.length = 0; // Add all bodies from the new world

            for (let i = 0; i < world.bodies.length; i++) {
              this.axisList.push(world.bodies[i]);
            } // Remove old handlers, if any

            world.removeEventListener("addBody", this._addBodyHandler);
            world.removeEventListener("removeBody", this._removeBodyHandler); // Add handlers to update the list of bodies.

            world.addEventListener("addBody", this._addBodyHandler);
            world.addEventListener("removeBody", this._removeBodyHandler);
            this.world = world;
            this.dirty = true;
          }
          /**
           * Collect all collision pairs
           */

          collisionPairs(world, p1, p2) {
            const bodies = this.axisList;
            const N = bodies.length;
            const axisIndex = this.axisIndex;
            let i;
            let j;

            if (this.dirty) {
              this.sortList();
              this.dirty = false;
            } // Look through the list

            for (i = 0; i !== N; i++) {
              const bi = bodies[i];

              for (j = i + 1; j < N; j++) {
                const bj = bodies[j];

                if (!this.needBroadphaseCollision(bi, bj)) {
                  continue;
                }

                if (!SAPBroadphase.checkBounds(bi, bj, axisIndex)) {
                  break;
                }

                this.intersectionTest(bi, bj, p1, p2);
              }
            }
          }

          sortList() {
            const axisList = this.axisList;
            const axisIndex = this.axisIndex;
            const N = axisList.length; // Update AABBs

            for (let i = 0; i !== N; i++) {
              const bi = axisList[i];

              if (bi.aabbNeedsUpdate) {
                bi.updateAABB();
              }
            } // Sort the list

            if (axisIndex === 0) {
              SAPBroadphase.insertionSortX(axisList);
            } else if (axisIndex === 1) {
              SAPBroadphase.insertionSortY(axisList);
            } else if (axisIndex === 2) {
              SAPBroadphase.insertionSortZ(axisList);
            }
          }
          /**
           * Computes the variance of the body positions and estimates the best axis to use.
           * Will automatically set property `axisIndex`.
           */

          autoDetectAxis() {
            let sumX = 0;
            let sumX2 = 0;
            let sumY = 0;
            let sumY2 = 0;
            let sumZ = 0;
            let sumZ2 = 0;
            const bodies = this.axisList;
            const N = bodies.length;
            const invN = 1 / N;

            for (let i = 0; i !== N; i++) {
              const b = bodies[i];
              const centerX = b.position.x;
              sumX += centerX;
              sumX2 += centerX * centerX;
              const centerY = b.position.y;
              sumY += centerY;
              sumY2 += centerY * centerY;
              const centerZ = b.position.z;
              sumZ += centerZ;
              sumZ2 += centerZ * centerZ;
            }

            const varianceX = sumX2 - sumX * sumX * invN;
            const varianceY = sumY2 - sumY * sumY * invN;
            const varianceZ = sumZ2 - sumZ * sumZ * invN;

            if (varianceX > varianceY) {
              if (varianceX > varianceZ) {
                this.axisIndex = 0;
              } else {
                this.axisIndex = 2;
              }
            } else if (varianceY > varianceZ) {
              this.axisIndex = 1;
            } else {
              this.axisIndex = 2;
            }
          }
          /**
           * Returns all the bodies within an AABB.
           * @param result An array to store resulting bodies in.
           */

          aabbQuery(world, aabb, result = []) {
            if (this.dirty) {
              this.sortList();
              this.dirty = false;
            }

            const axisIndex = this.axisIndex;
            let axis = "x";

            if (axisIndex === 1) {
              axis = "y";
            }

            if (axisIndex === 2) {
              axis = "z";
            }

            const axisList = this.axisList;
            aabb.lowerBound[axis];
            aabb.upperBound[axis];

            for (let i = 0; i < axisList.length; i++) {
              const b = axisList[i];

              if (b.aabbNeedsUpdate) {
                b.updateAABB();
              }

              if (b.aabb.overlaps(aabb)) {
                result.push(b);
              }
            }

            return result;
          }
        }

        class Utils {
          /**
           * Extend an options object with default values.
           * @param options The options object. May be falsy: in this case, a new object is created and returned.
           * @param defaults An object containing default values.
           * @return The modified options object.
           */
          static defaults(options = {}, defaults) {
            for (let key in defaults) {
              if (!(key in options)) {
                options[key] = defaults[key];
              }
            }

            return options;
          }
        }

        /**
         * Constraint base class
         */
        class Constraint {
          /**
           * Equations to be solved in this constraint.
           */

          /**
           * Body A.
           */

          /**
           * Body B.
           */

          /**
           * Set to false if you don't want the bodies to collide when they are connected.
           */
          constructor(bodyA, bodyB, options = {}) {
            this.equations = void 0;
            this.bodyA = void 0;
            this.bodyB = void 0;
            this.id = void 0;
            this.collideConnected = void 0;
            options = Utils.defaults(options, {
              collideConnected: true,
              wakeUpBodies: true,
            });
            this.equations = [];
            this.bodyA = bodyA;
            this.bodyB = bodyB;
            this.id = Constraint.idCounter++;
            this.collideConnected = options.collideConnected;

            if (options.wakeUpBodies) {
              if (bodyA) {
                bodyA.wakeUp();
              }

              if (bodyB) {
                bodyB.wakeUp();
              }
            }
          }
          /**
           * Update all the equations with data.
           */

          update() {
            throw new Error(
              "method update() not implmemented in this Constraint subclass!"
            );
          }
          /**
           * Enables all equations in the constraint.
           */

          enable() {
            const eqs = this.equations;

            for (let i = 0; i < eqs.length; i++) {
              eqs[i].enabled = true;
            }
          }
          /**
           * Disables all equations in the constraint.
           */

          disable() {
            const eqs = this.equations;

            for (let i = 0; i < eqs.length; i++) {
              eqs[i].enabled = false;
            }
          }
        }
        Constraint.idCounter = 0;

        /**
         * An element containing 6 entries, 3 spatial and 3 rotational degrees of freedom.
         */

        class JacobianElement {
          /**
           * spatial
           */

          /**
           * rotational
           */
          constructor() {
            this.spatial = void 0;
            this.rotational = void 0;
            this.spatial = new Vec3();
            this.rotational = new Vec3();
          }
          /**
           * Multiply with other JacobianElement
           */

          multiplyElement(element) {
            return (
              element.spatial.dot(this.spatial) +
              element.rotational.dot(this.rotational)
            );
          }
          /**
           * Multiply with two vectors
           */

          multiplyVectors(spatial, rotational) {
            return spatial.dot(this.spatial) + rotational.dot(this.rotational);
          }
        }

        /**
         * Equation base class.
         *
         * `a`, `b` and `eps` are {@link https://www8.cs.umu.se/kurser/5DV058/VT15/lectures/SPOOKlabnotes.pdf SPOOK} parameters that default to `0.0`. See {@link https://github.com/schteppe/cannon.js/issues/238#issuecomment-147172327 this exchange} for more details on Cannon's physics implementation.
         */
        class Equation {
          /**
           * Minimum (read: negative max) force to be applied by the constraint.
           */

          /**
           * Maximum (read: positive max) force to be applied by the constraint.
           */

          /**
           * SPOOK parameter
           */

          /**
           * SPOOK parameter
           */

          /**
           * SPOOK parameter
           */

          /**
           * A number, proportional to the force added to the bodies.
           */
          constructor(bi, bj, minForce = -1e6, maxForce = 1e6) {
            this.id = void 0;
            this.minForce = void 0;
            this.maxForce = void 0;
            this.bi = void 0;
            this.bj = void 0;
            this.si = void 0;
            this.sj = void 0;
            this.a = void 0;
            this.b = void 0;
            this.eps = void 0;
            this.jacobianElementA = void 0;
            this.jacobianElementB = void 0;
            this.enabled = void 0;
            this.multiplier = void 0;
            this.id = Equation.idCounter++;
            this.minForce = minForce;
            this.maxForce = maxForce;
            this.bi = bi;
            this.bj = bj;
            this.a = 0.0; // SPOOK parameter

            this.b = 0.0; // SPOOK parameter

            this.eps = 0.0; // SPOOK parameter

            this.jacobianElementA = new JacobianElement();
            this.jacobianElementB = new JacobianElement();
            this.enabled = true;
            this.multiplier = 0;
            this.setSpookParams(1e7, 4, 1 / 60); // Set typical spook params
          }
          /**
           * Recalculates a, b, and eps.
           *
           * The Equation constructor sets typical SPOOK parameters as such:
           * * `stiffness` = 1e7
           * * `relaxation` = 4
           * * `timeStep`= 1 / 60, _note the hardcoded refresh rate._
           */

          setSpookParams(stiffness, relaxation, timeStep) {
            const d = relaxation;
            const k = stiffness;
            const h = timeStep;
            this.a = 4.0 / (h * (1 + 4 * d));
            this.b = (4.0 * d) / (1 + 4 * d);
            this.eps = 4.0 / (h * h * k * (1 + 4 * d));
          }
          /**
           * Computes the right hand side of the SPOOK equation
           */

          computeB(a, b, h) {
            const GW = this.computeGW();
            const Gq = this.computeGq();
            const GiMf = this.computeGiMf();
            return -Gq * a - GW * b - GiMf * h;
          }
          /**
           * Computes G*q, where q are the generalized body coordinates
           */

          computeGq() {
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            const bi = this.bi;
            const bj = this.bj;
            const xi = bi.position;
            const xj = bj.position;
            return GA.spatial.dot(xi) + GB.spatial.dot(xj);
          }
          /**
           * Computes G*W, where W are the body velocities
           */

          computeGW() {
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            const bi = this.bi;
            const bj = this.bj;
            const vi = bi.velocity;
            const vj = bj.velocity;
            const wi = bi.angularVelocity;
            const wj = bj.angularVelocity;
            return GA.multiplyVectors(vi, wi) + GB.multiplyVectors(vj, wj);
          }
          /**
           * Computes G*Wlambda, where W are the body velocities
           */

          computeGWlambda() {
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            const bi = this.bi;
            const bj = this.bj;
            const vi = bi.vlambda;
            const vj = bj.vlambda;
            const wi = bi.wlambda;
            const wj = bj.wlambda;
            return GA.multiplyVectors(vi, wi) + GB.multiplyVectors(vj, wj);
          }
          /**
           * Computes G*inv(M)*f, where M is the mass matrix with diagonal blocks for each body, and f are the forces on the bodies.
           */

          computeGiMf() {
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            const bi = this.bi;
            const bj = this.bj;
            const fi = bi.force;
            const ti = bi.torque;
            const fj = bj.force;
            const tj = bj.torque;
            const invMassi = bi.invMassSolve;
            const invMassj = bj.invMassSolve;
            fi.scale(invMassi, iMfi);
            fj.scale(invMassj, iMfj);
            bi.invInertiaWorldSolve.vmult(ti, invIi_vmult_taui);
            bj.invInertiaWorldSolve.vmult(tj, invIj_vmult_tauj);
            return (
              GA.multiplyVectors(iMfi, invIi_vmult_taui) +
              GB.multiplyVectors(iMfj, invIj_vmult_tauj)
            );
          }
          /**
           * Computes G*inv(M)*G'
           */

          computeGiMGt() {
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            const bi = this.bi;
            const bj = this.bj;
            const invMassi = bi.invMassSolve;
            const invMassj = bj.invMassSolve;
            const invIi = bi.invInertiaWorldSolve;
            const invIj = bj.invInertiaWorldSolve;
            let result = invMassi + invMassj;
            invIi.vmult(GA.rotational, tmp);
            result += tmp.dot(GA.rotational);
            invIj.vmult(GB.rotational, tmp);
            result += tmp.dot(GB.rotational);
            return result;
          }
          /**
           * Add constraint velocity to the bodies.
           */

          addToWlambda(deltalambda) {
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            const bi = this.bi;
            const bj = this.bj;
            const temp = addToWlambda_temp; // Add to linear velocity
            // v_lambda += inv(M) * delta_lamba * G

            bi.vlambda.addScaledVector(
              bi.invMassSolve * deltalambda,
              GA.spatial,
              bi.vlambda
            );
            bj.vlambda.addScaledVector(
              bj.invMassSolve * deltalambda,
              GB.spatial,
              bj.vlambda
            ); // Add to angular velocity

            bi.invInertiaWorldSolve.vmult(GA.rotational, temp);
            bi.wlambda.addScaledVector(deltalambda, temp, bi.wlambda);
            bj.invInertiaWorldSolve.vmult(GB.rotational, temp);
            bj.wlambda.addScaledVector(deltalambda, temp, bj.wlambda);
          }
          /**
           * Compute the denominator part of the SPOOK equation: C = G*inv(M)*G' + eps
           */

          computeC() {
            return this.computeGiMGt() + this.eps;
          }
        }
        Equation.idCounter = 0;
        const iMfi = new Vec3();
        const iMfj = new Vec3();
        const invIi_vmult_taui = new Vec3();
        const invIj_vmult_tauj = new Vec3();
        const tmp = new Vec3();
        const addToWlambda_temp = new Vec3();

        /**
         * Contact/non-penetration constraint equation
         */
        class ContactEquation extends Equation {
          /**
           * "bounciness": u1 = -e*u0
           */

          /**
           * World-oriented vector that goes from the center of bi to the contact point.
           */

          /**
           * World-oriented vector that starts in body j position and goes to the contact point.
           */

          /**
           * Contact normal, pointing out of body i.
           */
          constructor(bodyA, bodyB, maxForce = 1e6) {
            super(bodyA, bodyB, 0, maxForce);
            this.restitution = void 0;
            this.ri = void 0;
            this.rj = void 0;
            this.ni = void 0;
            this.restitution = 0.0;
            this.ri = new Vec3();
            this.rj = new Vec3();
            this.ni = new Vec3();
          }

          computeB(h) {
            const a = this.a;
            const b = this.b;
            const bi = this.bi;
            const bj = this.bj;
            const ri = this.ri;
            const rj = this.rj;
            const rixn = ContactEquation_computeB_temp1;
            const rjxn = ContactEquation_computeB_temp2;
            const vi = bi.velocity;
            const wi = bi.angularVelocity;
            bi.force;
            bi.torque;
            const vj = bj.velocity;
            const wj = bj.angularVelocity;
            bj.force;
            bj.torque;
            const penetrationVec = ContactEquation_computeB_temp3;
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            const n = this.ni; // Caluclate cross products

            ri.cross(n, rixn);
            rj.cross(n, rjxn); // g = xj+rj -(xi+ri)
            // G = [ -ni  -rixn  ni  rjxn ]

            n.negate(GA.spatial);
            rixn.negate(GA.rotational);
            GB.spatial.copy(n);
            GB.rotational.copy(rjxn); // Calculate the penetration vector

            penetrationVec.copy(bj.position);
            penetrationVec.vadd(rj, penetrationVec);
            penetrationVec.vsub(bi.position, penetrationVec);
            penetrationVec.vsub(ri, penetrationVec);
            const g = n.dot(penetrationVec); // Compute iteration

            const ePlusOne = this.restitution + 1;
            const GW =
              ePlusOne * vj.dot(n) -
              ePlusOne * vi.dot(n) +
              wj.dot(rjxn) -
              wi.dot(rixn);
            const GiMf = this.computeGiMf();
            const B = -g * a - GW * b - h * GiMf;
            return B;
          }
          /**
           * Get the current relative velocity in the contact point.
           */

          getImpactVelocityAlongNormal() {
            const vi = ContactEquation_getImpactVelocityAlongNormal_vi;
            const vj = ContactEquation_getImpactVelocityAlongNormal_vj;
            const xi = ContactEquation_getImpactVelocityAlongNormal_xi;
            const xj = ContactEquation_getImpactVelocityAlongNormal_xj;
            const relVel = ContactEquation_getImpactVelocityAlongNormal_relVel;
            this.bi.position.vadd(this.ri, xi);
            this.bj.position.vadd(this.rj, xj);
            this.bi.getVelocityAtWorldPoint(xi, vi);
            this.bj.getVelocityAtWorldPoint(xj, vj);
            vi.vsub(vj, relVel);
            return this.ni.dot(relVel);
          }
        }
        const ContactEquation_computeB_temp1 = new Vec3(); // Temp vectors

        const ContactEquation_computeB_temp2 = new Vec3();
        const ContactEquation_computeB_temp3 = new Vec3();
        const ContactEquation_getImpactVelocityAlongNormal_vi = new Vec3();
        const ContactEquation_getImpactVelocityAlongNormal_vj = new Vec3();
        const ContactEquation_getImpactVelocityAlongNormal_xi = new Vec3();
        const ContactEquation_getImpactVelocityAlongNormal_xj = new Vec3();
        const ContactEquation_getImpactVelocityAlongNormal_relVel = new Vec3();

        /**
         * Connects two bodies at given offset points.
         * @example
         *     const bodyA = new Body({ mass: 1 })
         *     const bodyB = new Body({ mass: 1 })
         *     bodyA.position.set(-1, 0, 0)
         *     bodyB.position.set(1, 0, 0)
         *     bodyA.addShape(shapeA)
         *     bodyB.addShape(shapeB)
         *     world.addBody(bodyA)
         *     world.addBody(bodyB)
         *     const localPivotA = new Vec3(1, 0, 0)
         *     const localPivotB = new Vec3(-1, 0, 0)
         *     const constraint = new PointToPointConstraint(bodyA, localPivotA, bodyB, localPivotB)
         *     world.addConstraint(constraint)
         */
        class PointToPointConstraint extends Constraint {
          /**
           * Pivot, defined locally in bodyA.
           */

          /**
           * Pivot, defined locally in bodyB.
           */

          /**
           * @param pivotA The point relative to the center of mass of bodyA which bodyA is constrained to.
           * @param bodyB Body that will be constrained in a similar way to the same point as bodyA. We will therefore get a link between bodyA and bodyB. If not specified, bodyA will be constrained to a static point.
           * @param pivotB The point relative to the center of mass of bodyB which bodyB is constrained to.
           * @param maxForce The maximum force that should be applied to constrain the bodies.
           */
          constructor(
            bodyA,
            pivotA = new Vec3(),
            bodyB,
            pivotB = new Vec3(),
            maxForce = 1e6
          ) {
            super(bodyA, bodyB);
            this.pivotA = void 0;
            this.pivotB = void 0;
            this.equationX = void 0;
            this.equationY = void 0;
            this.equationZ = void 0;
            this.pivotA = pivotA.clone();
            this.pivotB = pivotB.clone();
            const x = (this.equationX = new ContactEquation(bodyA, bodyB));
            const y = (this.equationY = new ContactEquation(bodyA, bodyB));
            const z = (this.equationZ = new ContactEquation(bodyA, bodyB)); // Equations to be fed to the solver

            this.equations.push(x, y, z); // Make the equations bidirectional

            x.minForce = y.minForce = z.minForce = -maxForce;
            x.maxForce = y.maxForce = z.maxForce = maxForce;
            x.ni.set(1, 0, 0);
            y.ni.set(0, 1, 0);
            z.ni.set(0, 0, 1);
          }

          update() {
            const bodyA = this.bodyA;
            const bodyB = this.bodyB;
            const x = this.equationX;
            const y = this.equationY;
            const z = this.equationZ; // Rotate the pivots to world space

            bodyA.quaternion.vmult(this.pivotA, x.ri);
            bodyB.quaternion.vmult(this.pivotB, x.rj);
            y.ri.copy(x.ri);
            y.rj.copy(x.rj);
            z.ri.copy(x.ri);
            z.rj.copy(x.rj);
          }
        }

        /**
         * Cone equation. Works to keep the given body world vectors aligned, or tilted within a given angle from each other.
         */
        class ConeEquation extends Equation {
          /**
           * Local axis in A
           */

          /**
           * Local axis in B
           */

          /**
           * The "cone angle" to keep
           */
          constructor(bodyA, bodyB, options = {}) {
            const maxForce =
              typeof options.maxForce !== "undefined" ? options.maxForce : 1e6;
            super(bodyA, bodyB, -maxForce, maxForce);
            this.axisA = void 0;
            this.axisB = void 0;
            this.angle = void 0;
            this.axisA = options.axisA
              ? options.axisA.clone()
              : new Vec3(1, 0, 0);
            this.axisB = options.axisB
              ? options.axisB.clone()
              : new Vec3(0, 1, 0);
            this.angle =
              typeof options.angle !== "undefined" ? options.angle : 0;
          }

          computeB(h) {
            const a = this.a;
            const b = this.b;
            const ni = this.axisA;
            const nj = this.axisB;
            const nixnj = tmpVec1$2;
            const njxni = tmpVec2$2;
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB; // Caluclate cross products

            ni.cross(nj, nixnj);
            nj.cross(ni, njxni); // The angle between two vector is:
            // cos(theta) = a * b / (length(a) * length(b) = { len(a) = len(b) = 1 } = a * b
            // g = a * b
            // gdot = (b x a) * wi + (a x b) * wj
            // G = [0 bxa 0 axb]
            // W = [vi wi vj wj]

            GA.rotational.copy(njxni);
            GB.rotational.copy(nixnj);
            const g = Math.cos(this.angle) - ni.dot(nj);
            const GW = this.computeGW();
            const GiMf = this.computeGiMf();
            const B = -g * a - GW * b - h * GiMf;
            return B;
          }
        }
        const tmpVec1$2 = new Vec3();
        const tmpVec2$2 = new Vec3();

        /**
         * Rotational constraint. Works to keep the local vectors orthogonal to each other in world space.
         */
        class RotationalEquation extends Equation {
          /**
           * World oriented rotational axis.
           */

          /**
           * World oriented rotational axis.
           */

          /**
           * maxAngle
           */
          constructor(bodyA, bodyB, options = {}) {
            const maxForce =
              typeof options.maxForce !== "undefined" ? options.maxForce : 1e6;
            super(bodyA, bodyB, -maxForce, maxForce);
            this.axisA = void 0;
            this.axisB = void 0;
            this.maxAngle = void 0;
            this.axisA = options.axisA
              ? options.axisA.clone()
              : new Vec3(1, 0, 0);
            this.axisB = options.axisB
              ? options.axisB.clone()
              : new Vec3(0, 1, 0);
            this.maxAngle = Math.PI / 2;
          }

          computeB(h) {
            const a = this.a;
            const b = this.b;
            const ni = this.axisA;
            const nj = this.axisB;
            const nixnj = tmpVec1$1;
            const njxni = tmpVec2$1;
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB; // Caluclate cross products

            ni.cross(nj, nixnj);
            nj.cross(ni, njxni); // g = ni * nj
            // gdot = (nj x ni) * wi + (ni x nj) * wj
            // G = [0 njxni 0 nixnj]
            // W = [vi wi vj wj]

            GA.rotational.copy(njxni);
            GB.rotational.copy(nixnj);
            const g = Math.cos(this.maxAngle) - ni.dot(nj);
            const GW = this.computeGW();
            const GiMf = this.computeGiMf();
            const B = -g * a - GW * b - h * GiMf;
            return B;
          }
        }
        const tmpVec1$1 = new Vec3();
        const tmpVec2$1 = new Vec3();

        /**
         * A Cone Twist constraint, useful for ragdolls.
         */
        class ConeTwistConstraint extends PointToPointConstraint {
          /**
           * The axis direction for the constraint of the body A.
           */

          /**
           * The axis direction for the constraint of the body B.
           */

          /**
           * The aperture angle of the cone.
           */

          /**
           * The twist angle of the joint.
           */
          constructor(bodyA, bodyB, options = {}) {
            const maxForce =
              typeof options.maxForce !== "undefined" ? options.maxForce : 1e6; // Set pivot point in between

            const pivotA = options.pivotA ? options.pivotA.clone() : new Vec3();
            const pivotB = options.pivotB ? options.pivotB.clone() : new Vec3();
            super(bodyA, pivotA, bodyB, pivotB, maxForce);
            this.axisA = void 0;
            this.axisB = void 0;
            this.angle = void 0;
            this.twistAngle = void 0;
            this.coneEquation = void 0;
            this.twistEquation = void 0;
            this.axisA = options.axisA ? options.axisA.clone() : new Vec3();
            this.axisB = options.axisB ? options.axisB.clone() : new Vec3();
            this.collideConnected = !!options.collideConnected;
            this.angle =
              typeof options.angle !== "undefined" ? options.angle : 0;
            const c = (this.coneEquation = new ConeEquation(
              bodyA,
              bodyB,
              options
            ));
            const t = (this.twistEquation = new RotationalEquation(
              bodyA,
              bodyB,
              options
            ));
            this.twistAngle =
              typeof options.twistAngle !== "undefined"
                ? options.twistAngle
                : 0; // Make the cone equation push the bodies toward the cone axis, not outward

            c.maxForce = 0;
            c.minForce = -maxForce; // Make the twist equation add torque toward the initial position

            t.maxForce = 0;
            t.minForce = -maxForce;
            this.equations.push(c, t);
          }

          update() {
            const bodyA = this.bodyA;
            const bodyB = this.bodyB;
            const cone = this.coneEquation;
            const twist = this.twistEquation;
            super.update(); // Update the axes to the cone constraint

            bodyA.vectorToWorldFrame(this.axisA, cone.axisA);
            bodyB.vectorToWorldFrame(this.axisB, cone.axisB); // Update the world axes in the twist constraint

            this.axisA.tangents(twist.axisA, twist.axisA);
            bodyA.vectorToWorldFrame(twist.axisA, twist.axisA);
            this.axisB.tangents(twist.axisB, twist.axisB);
            bodyB.vectorToWorldFrame(twist.axisB, twist.axisB);
            cone.angle = this.angle;
            twist.maxAngle = this.twistAngle;
          }
        }

        /**
         * Constrains two bodies to be at a constant distance from each others center of mass.
         */
        class DistanceConstraint extends Constraint {
          /**
           * The distance to keep. If undefined, it will be set to the current distance between bodyA and bodyB
           */

          /**
           * @param distance The distance to keep. If undefined, it will be set to the current distance between bodyA and bodyB.
           * @param maxForce The maximum force that should be applied to constrain the bodies.
           */
          constructor(bodyA, bodyB, distance, maxForce = 1e6) {
            super(bodyA, bodyB);
            this.distance = void 0;
            this.distanceEquation = void 0;

            if (typeof distance === "undefined") {
              distance = bodyA.position.distanceTo(bodyB.position);
            }

            this.distance = distance;
            const eq = (this.distanceEquation = new ContactEquation(
              bodyA,
              bodyB
            ));
            this.equations.push(eq); // Make it bidirectional

            eq.minForce = -maxForce;
            eq.maxForce = maxForce;
          }
          /**
           * update
           */

          update() {
            const bodyA = this.bodyA;
            const bodyB = this.bodyB;
            const eq = this.distanceEquation;
            const halfDist = this.distance * 0.5;
            const normal = eq.ni;
            bodyB.position.vsub(bodyA.position, normal);
            normal.normalize();
            normal.scale(halfDist, eq.ri);
            normal.scale(-halfDist, eq.rj);
          }
        }

        /**
         * Lock constraint. Will remove all degrees of freedom between the bodies.
         */
        class LockConstraint extends PointToPointConstraint {
          constructor(bodyA, bodyB, options = {}) {
            const maxForce =
              typeof options.maxForce !== "undefined" ? options.maxForce : 1e6; // Set pivot point in between

            const pivotA = new Vec3();
            const pivotB = new Vec3();
            const halfWay = new Vec3();
            bodyA.position.vadd(bodyB.position, halfWay);
            halfWay.scale(0.5, halfWay);
            bodyB.pointToLocalFrame(halfWay, pivotB);
            bodyA.pointToLocalFrame(halfWay, pivotA); // The point-to-point constraint will keep a point shared between the bodies

            super(bodyA, pivotA, bodyB, pivotB, maxForce); // Store initial rotation of the bodies as unit vectors in the local body spaces

            this.xA = void 0;
            this.xB = void 0;
            this.yA = void 0;
            this.yB = void 0;
            this.zA = void 0;
            this.zB = void 0;
            this.rotationalEquation1 = void 0;
            this.rotationalEquation2 = void 0;
            this.rotationalEquation3 = void 0;
            this.motorEquation = void 0;
            this.xA = bodyA.vectorToLocalFrame(Vec3.UNIT_X);
            this.xB = bodyB.vectorToLocalFrame(Vec3.UNIT_X);
            this.yA = bodyA.vectorToLocalFrame(Vec3.UNIT_Y);
            this.yB = bodyB.vectorToLocalFrame(Vec3.UNIT_Y);
            this.zA = bodyA.vectorToLocalFrame(Vec3.UNIT_Z);
            this.zB = bodyB.vectorToLocalFrame(Vec3.UNIT_Z); // ...and the following rotational equations will keep all rotational DOF's in place

            const r1 = (this.rotationalEquation1 = new RotationalEquation(
              bodyA,
              bodyB,
              options
            ));
            const r2 = (this.rotationalEquation2 = new RotationalEquation(
              bodyA,
              bodyB,
              options
            ));
            const r3 = (this.rotationalEquation3 = new RotationalEquation(
              bodyA,
              bodyB,
              options
            ));
            this.equations.push(r1, r2, r3);
          }
          /**
           * update
           */

          update() {
            const bodyA = this.bodyA;
            const bodyB = this.bodyB;
            this.motorEquation;
            const r1 = this.rotationalEquation1;
            const r2 = this.rotationalEquation2;
            const r3 = this.rotationalEquation3;
            super.update(); // These vector pairs must be orthogonal

            bodyA.vectorToWorldFrame(this.xA, r1.axisA);
            bodyB.vectorToWorldFrame(this.yB, r1.axisB);
            bodyA.vectorToWorldFrame(this.yA, r2.axisA);
            bodyB.vectorToWorldFrame(this.zB, r2.axisB);
            bodyA.vectorToWorldFrame(this.zA, r3.axisA);
            bodyB.vectorToWorldFrame(this.xB, r3.axisB);
          }
        }

        /**
         * Rotational motor constraint. Tries to keep the relative angular velocity of the bodies to a given value.
         */
        class RotationalMotorEquation extends Equation {
          /**
           * World oriented rotational axis.
           */

          /**
           * World oriented rotational axis.
           */

          /**
           * Motor velocity.
           */
          constructor(bodyA, bodyB, maxForce = 1e6) {
            super(bodyA, bodyB, -maxForce, maxForce);
            this.axisA = void 0;
            this.axisB = void 0;
            this.targetVelocity = void 0;
            this.axisA = new Vec3();
            this.axisB = new Vec3();
            this.targetVelocity = 0;
          }

          computeB(h) {
            this.a;
            const b = this.b;
            this.bi;
            this.bj;
            const axisA = this.axisA;
            const axisB = this.axisB;
            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB; // g = 0
            // gdot = axisA * wi - axisB * wj
            // gdot = G * W = G * [vi wi vj wj]
            // =>
            // G = [0 axisA 0 -axisB]

            GA.rotational.copy(axisA);
            axisB.negate(GB.rotational);
            const GW = this.computeGW() - this.targetVelocity;
            const GiMf = this.computeGiMf();
            const B = -GW * b - h * GiMf;
            return B;
          }
        }

        /**
         * Hinge constraint. Think of it as a door hinge. It tries to keep the door in the correct place and with the correct orientation.
         */
        class HingeConstraint extends PointToPointConstraint {
          /**
           * Rotation axis, defined locally in bodyA.
           */

          /**
           * Rotation axis, defined locally in bodyB.
           */
          constructor(bodyA, bodyB, options = {}) {
            const maxForce =
              typeof options.maxForce !== "undefined" ? options.maxForce : 1e6;
            const pivotA = options.pivotA ? options.pivotA.clone() : new Vec3();
            const pivotB = options.pivotB ? options.pivotB.clone() : new Vec3();
            super(bodyA, pivotA, bodyB, pivotB, maxForce);
            this.axisA = void 0;
            this.axisB = void 0;
            this.rotationalEquation1 = void 0;
            this.rotationalEquation2 = void 0;
            this.motorEquation = void 0;
            const axisA = (this.axisA = options.axisA
              ? options.axisA.clone()
              : new Vec3(1, 0, 0));
            axisA.normalize();
            const axisB = (this.axisB = options.axisB
              ? options.axisB.clone()
              : new Vec3(1, 0, 0));
            axisB.normalize();
            this.collideConnected = !!options.collideConnected;
            const rotational1 = (this.rotationalEquation1 =
              new RotationalEquation(bodyA, bodyB, options));
            const rotational2 = (this.rotationalEquation2 =
              new RotationalEquation(bodyA, bodyB, options));
            const motor = (this.motorEquation = new RotationalMotorEquation(
              bodyA,
              bodyB,
              maxForce
            ));
            motor.enabled = false; // Not enabled by default
            // Equations to be fed to the solver

            this.equations.push(rotational1, rotational2, motor);
          }
          /**
           * enableMotor
           */

          enableMotor() {
            this.motorEquation.enabled = true;
          }
          /**
           * disableMotor
           */

          disableMotor() {
            this.motorEquation.enabled = false;
          }
          /**
           * setMotorSpeed
           */

          setMotorSpeed(speed) {
            this.motorEquation.targetVelocity = speed;
          }
          /**
           * setMotorMaxForce
           */

          setMotorMaxForce(maxForce) {
            this.motorEquation.maxForce = maxForce;
            this.motorEquation.minForce = -maxForce;
          }
          /**
           * update
           */

          update() {
            const bodyA = this.bodyA;
            const bodyB = this.bodyB;
            const motor = this.motorEquation;
            const r1 = this.rotationalEquation1;
            const r2 = this.rotationalEquation2;
            const worldAxisA = HingeConstraint_update_tmpVec1;
            const worldAxisB = HingeConstraint_update_tmpVec2;
            const axisA = this.axisA;
            const axisB = this.axisB;
            super.update(); // Get world axes

            bodyA.quaternion.vmult(axisA, worldAxisA);
            bodyB.quaternion.vmult(axisB, worldAxisB);
            worldAxisA.tangents(r1.axisA, r2.axisA);
            r1.axisB.copy(worldAxisB);
            r2.axisB.copy(worldAxisB);

            if (this.motorEquation.enabled) {
              bodyA.quaternion.vmult(this.axisA, motor.axisA);
              bodyB.quaternion.vmult(this.axisB, motor.axisB);
            }
          }
        }
        const HingeConstraint_update_tmpVec1 = new Vec3();
        const HingeConstraint_update_tmpVec2 = new Vec3();

        /**
         * Constrains the slipping in a contact along a tangent
         */
        class FrictionEquation extends Equation {
          // Tangent

          /**
           * @param slipForce should be +-F_friction = +-mu * F_normal = +-mu * m * g
           */
          constructor(bodyA, bodyB, slipForce) {
            super(bodyA, bodyB, -slipForce, slipForce);
            this.ri = void 0;
            this.rj = void 0;
            this.t = void 0;
            this.ri = new Vec3();
            this.rj = new Vec3();
            this.t = new Vec3();
          }

          computeB(h) {
            this.a;
            const b = this.b;
            this.bi;
            this.bj;
            const ri = this.ri;
            const rj = this.rj;
            const rixt = FrictionEquation_computeB_temp1;
            const rjxt = FrictionEquation_computeB_temp2;
            const t = this.t; // Caluclate cross products

            ri.cross(t, rixt);
            rj.cross(t, rjxt); // G = [-t -rixt t rjxt]
            // And remember, this is a pure velocity constraint, g is always zero!

            const GA = this.jacobianElementA;
            const GB = this.jacobianElementB;
            t.negate(GA.spatial);
            rixt.negate(GA.rotational);
            GB.spatial.copy(t);
            GB.rotational.copy(rjxt);
            const GW = this.computeGW();
            const GiMf = this.computeGiMf();
            const B = -GW * b - h * GiMf;
            return B;
          }
        }
        const FrictionEquation_computeB_temp1 = new Vec3();
        const FrictionEquation_computeB_temp2 = new Vec3();

        /**
         * Defines what happens when two materials meet.
         * @todo Refactor materials to materialA and materialB
         */
        class ContactMaterial {
          /**
           * Identifier of this material.
           */

          /**
           * Participating materials.
           */

          /**
           * Friction coefficient.
           * @default 0.3
           */

          /**
           * Restitution coefficient.
           * @default 0.3
           */

          /**
           * Stiffness of the produced contact equations.
           * @default 1e7
           */

          /**
           * Relaxation time of the produced contact equations.
           * @default 3
           */

          /**
           * Stiffness of the produced friction equations.
           * @default 1e7
           */

          /**
           * Relaxation time of the produced friction equations
           * @default 3
           */
          constructor(m1, m2, options) {
            this.id = void 0;
            this.materials = void 0;
            this.friction = void 0;
            this.restitution = void 0;
            this.contactEquationStiffness = void 0;
            this.contactEquationRelaxation = void 0;
            this.frictionEquationStiffness = void 0;
            this.frictionEquationRelaxation = void 0;
            options = Utils.defaults(options, {
              friction: 0.3,
              restitution: 0.3,
              contactEquationStiffness: 1e7,
              contactEquationRelaxation: 3,
              frictionEquationStiffness: 1e7,
              frictionEquationRelaxation: 3,
            });
            this.id = ContactMaterial.idCounter++;
            this.materials = [m1, m2];
            this.friction = options.friction;
            this.restitution = options.restitution;
            this.contactEquationStiffness = options.contactEquationStiffness;
            this.contactEquationRelaxation = options.contactEquationRelaxation;
            this.frictionEquationStiffness = options.frictionEquationStiffness;
            this.frictionEquationRelaxation =
              options.frictionEquationRelaxation;
          }
        }
        ContactMaterial.idCounter = 0;

        /**
         * Defines a physics material.
         */
        class Material {
          /**
           * Material name.
           * If options is a string, name will be set to that string.
           * @todo Deprecate this
           */

          /** Material id. */

          /**
           * Friction for this material.
           * If non-negative, it will be used instead of the friction given by ContactMaterials. If there's no matching ContactMaterial, the value from `defaultContactMaterial` in the World will be used.
           */

          /**
           * Restitution for this material.
           * If non-negative, it will be used instead of the restitution given by ContactMaterials. If there's no matching ContactMaterial, the value from `defaultContactMaterial` in the World will be used.
           */
          constructor(options = {}) {
            this.name = void 0;
            this.id = void 0;
            this.friction = void 0;
            this.restitution = void 0;
            let name = ""; // Backwards compatibility fix

            if (typeof options === "string") {
              //console.warn(`Passing a string to MaterialOptions is deprecated, and has no effect`)
              name = options;
              options = {};
            }

            this.name = name;
            this.id = Material.idCounter++;
            this.friction =
              typeof options.friction !== "undefined" ? options.friction : -1;
            this.restitution =
              typeof options.restitution !== "undefined"
                ? options.restitution
                : -1;
          }
        }
        Material.idCounter = 0;

        /**
         * A spring, connecting two bodies.
         * @example
         *     const spring = new Spring(boxBody, sphereBody, {
         *       restLength: 0,
         *       stiffness: 50,
         *       damping: 1,
         *     })
         *
         *     // Compute the force after each step
         *     world.addEventListener('postStep', (event) => {
         *       spring.applyForce()
         *     })
         */
        class Spring {
          /**
           * Rest length of the spring. A number > 0.
           * @default 1
           */

          /**
           * Stiffness of the spring. A number >= 0.
           * @default 100
           */

          /**
           * Damping of the spring. A number >= 0.
           * @default 1
           */

          /**
           * First connected body.
           */

          /**
           * Second connected body.
           */

          /**
           * Anchor for bodyA in local bodyA coordinates.
           * Where to hook the spring to body A, in local body coordinates.
           * @default new Vec3()
           */

          /**
           * Anchor for bodyB in local bodyB coordinates.
           * Where to hook the spring to body B, in local body coordinates.
           * @default new Vec3()
           */
          constructor(bodyA, bodyB, options = {}) {
            this.restLength = void 0;
            this.stiffness = void 0;
            this.damping = void 0;
            this.bodyA = void 0;
            this.bodyB = void 0;
            this.localAnchorA = void 0;
            this.localAnchorB = void 0;
            this.restLength =
              typeof options.restLength === "number" ? options.restLength : 1;
            this.stiffness = options.stiffness || 100;
            this.damping = options.damping || 1;
            this.bodyA = bodyA;
            this.bodyB = bodyB;
            this.localAnchorA = new Vec3();
            this.localAnchorB = new Vec3();

            if (options.localAnchorA) {
              this.localAnchorA.copy(options.localAnchorA);
            }

            if (options.localAnchorB) {
              this.localAnchorB.copy(options.localAnchorB);
            }

            if (options.worldAnchorA) {
              this.setWorldAnchorA(options.worldAnchorA);
            }

            if (options.worldAnchorB) {
              this.setWorldAnchorB(options.worldAnchorB);
            }
          }
          /**
           * Set the anchor point on body A, using world coordinates.
           */

          setWorldAnchorA(worldAnchorA) {
            this.bodyA.pointToLocalFrame(worldAnchorA, this.localAnchorA);
          }
          /**
           * Set the anchor point on body B, using world coordinates.
           */

          setWorldAnchorB(worldAnchorB) {
            this.bodyB.pointToLocalFrame(worldAnchorB, this.localAnchorB);
          }
          /**
           * Get the anchor point on body A, in world coordinates.
           * @param result The vector to store the result in.
           */

          getWorldAnchorA(result) {
            this.bodyA.pointToWorldFrame(this.localAnchorA, result);
          }
          /**
           * Get the anchor point on body B, in world coordinates.
           * @param result The vector to store the result in.
           */

          getWorldAnchorB(result) {
            this.bodyB.pointToWorldFrame(this.localAnchorB, result);
          }
          /**
           * Apply the spring force to the connected bodies.
           */

          applyForce() {
            const k = this.stiffness;
            const d = this.damping;
            const l = this.restLength;
            const bodyA = this.bodyA;
            const bodyB = this.bodyB;
            const r = applyForce_r;
            const r_unit = applyForce_r_unit;
            const u = applyForce_u;
            const f = applyForce_f;
            const tmp = applyForce_tmp;
            const worldAnchorA = applyForce_worldAnchorA;
            const worldAnchorB = applyForce_worldAnchorB;
            const ri = applyForce_ri;
            const rj = applyForce_rj;
            const ri_x_f = applyForce_ri_x_f;
            const rj_x_f = applyForce_rj_x_f; // Get world anchors

            this.getWorldAnchorA(worldAnchorA);
            this.getWorldAnchorB(worldAnchorB); // Get offset points

            worldAnchorA.vsub(bodyA.position, ri);
            worldAnchorB.vsub(bodyB.position, rj); // Compute distance vector between world anchor points

            worldAnchorB.vsub(worldAnchorA, r);
            const rlen = r.length();
            r_unit.copy(r);
            r_unit.normalize(); // Compute relative velocity of the anchor points, u

            bodyB.velocity.vsub(bodyA.velocity, u); // Add rotational velocity

            bodyB.angularVelocity.cross(rj, tmp);
            u.vadd(tmp, u);
            bodyA.angularVelocity.cross(ri, tmp);
            u.vsub(tmp, u); // F = - k * ( x - L ) - D * ( u )

            r_unit.scale(-k * (rlen - l) - d * u.dot(r_unit), f); // Add forces to bodies

            bodyA.force.vsub(f, bodyA.force);
            bodyB.force.vadd(f, bodyB.force); // Angular force

            ri.cross(f, ri_x_f);
            rj.cross(f, rj_x_f);
            bodyA.torque.vsub(ri_x_f, bodyA.torque);
            bodyB.torque.vadd(rj_x_f, bodyB.torque);
          }
        }
        const applyForce_r = new Vec3();
        const applyForce_r_unit = new Vec3();
        const applyForce_u = new Vec3();
        const applyForce_f = new Vec3();
        const applyForce_worldAnchorA = new Vec3();
        const applyForce_worldAnchorB = new Vec3();
        const applyForce_ri = new Vec3();
        const applyForce_rj = new Vec3();
        const applyForce_ri_x_f = new Vec3();
        const applyForce_rj_x_f = new Vec3();
        const applyForce_tmp = new Vec3();

        /**
         * WheelInfo
         */
        class WheelInfo {
          /**
           * Max travel distance of the suspension, in meters.
           * @default 1
           */

          /**
           * Speed to apply to the wheel rotation when the wheel is sliding.
           * @default -0.1
           */

          /**
           * If the customSlidingRotationalSpeed should be used.
           * @default false
           */

          /**
           * sliding
           */

          /**
           * Connection point, defined locally in the chassis body frame.
           */

          /**
           * chassisConnectionPointWorld
           */

          /**
           * directionLocal
           */

          /**
           * directionWorld
           */

          /**
           * axleLocal
           */

          /**
           * axleWorld
           */

          /**
           * suspensionRestLength
           * @default 1
           */

          /**
           * suspensionMaxLength
           * @default 2
           */

          /**
           * radius
           * @default 1
           */

          /**
           * suspensionStiffness
           * @default 100
           */

          /**
           * dampingCompression
           * @default 10
           */

          /**
           * dampingRelaxation
           * @default 10
           */

          /**
           * frictionSlip
           * @default 10.5
           */

          /** forwardAcceleration */

          /** sideAcceleration */

          /**
           * steering
           * @default 0
           */

          /**
           * Rotation value, in radians.
           * @default 0
           */

          /**
           * deltaRotation
           * @default 0
           */

          /**
           * rollInfluence
           * @default 0.01
           */

          /**
           * maxSuspensionForce
           */

          /**
           * engineForce
           */

          /**
           * brake
           */

          /**
           * isFrontWheel
           * @default true
           */

          /**
           * clippedInvContactDotSuspension
           * @default 1
           */

          /**
           * suspensionRelativeVelocity
           * @default 0
           */

          /**
           * suspensionForce
           * @default 0
           */

          /**
           * slipInfo
           */

          /**
           * skidInfo
           * @default 0
           */

          /**
           * suspensionLength
           * @default 0
           */

          /**
           * sideImpulse
           */

          /**
           * forwardImpulse
           */

          /**
           * The result from raycasting.
           */

          /**
           * Wheel world transform.
           */

          /**
           * isInContact
           */
          constructor(options = {}) {
            this.maxSuspensionTravel = void 0;
            this.customSlidingRotationalSpeed = void 0;
            this.useCustomSlidingRotationalSpeed = void 0;
            this.sliding = void 0;
            this.chassisConnectionPointLocal = void 0;
            this.chassisConnectionPointWorld = void 0;
            this.directionLocal = void 0;
            this.directionWorld = void 0;
            this.axleLocal = void 0;
            this.axleWorld = void 0;
            this.suspensionRestLength = void 0;
            this.suspensionMaxLength = void 0;
            this.radius = void 0;
            this.suspensionStiffness = void 0;
            this.dampingCompression = void 0;
            this.dampingRelaxation = void 0;
            this.frictionSlip = void 0;
            this.forwardAcceleration = void 0;
            this.sideAcceleration = void 0;
            this.steering = void 0;
            this.rotation = void 0;
            this.deltaRotation = void 0;
            this.rollInfluence = void 0;
            this.maxSuspensionForce = void 0;
            this.engineForce = void 0;
            this.brake = void 0;
            this.isFrontWheel = void 0;
            this.clippedInvContactDotSuspension = void 0;
            this.suspensionRelativeVelocity = void 0;
            this.suspensionForce = void 0;
            this.slipInfo = void 0;
            this.skidInfo = void 0;
            this.suspensionLength = void 0;
            this.sideImpulse = void 0;
            this.forwardImpulse = void 0;
            this.raycastResult = void 0;
            this.worldTransform = void 0;
            this.isInContact = void 0;
            options = Utils.defaults(options, {
              chassisConnectionPointLocal: new Vec3(),
              chassisConnectionPointWorld: new Vec3(),
              directionLocal: new Vec3(),
              directionWorld: new Vec3(),
              axleLocal: new Vec3(),
              axleWorld: new Vec3(),
              suspensionRestLength: 1,
              suspensionMaxLength: 2,
              radius: 1,
              suspensionStiffness: 100,
              dampingCompression: 10,
              dampingRelaxation: 10,
              frictionSlip: 10.5,
              forwardAcceleration: 1,
              sideAcceleration: 1,
              steering: 0,
              rotation: 0,
              deltaRotation: 0,
              rollInfluence: 0.01,
              maxSuspensionForce: Number.MAX_VALUE,
              isFrontWheel: true,
              clippedInvContactDotSuspension: 1,
              suspensionRelativeVelocity: 0,
              suspensionForce: 0,
              slipInfo: 0,
              skidInfo: 0,
              suspensionLength: 0,
              maxSuspensionTravel: 1,
              useCustomSlidingRotationalSpeed: false,
              customSlidingRotationalSpeed: -0.1,
            });
            this.maxSuspensionTravel = options.maxSuspensionTravel;
            this.customSlidingRotationalSpeed =
              options.customSlidingRotationalSpeed;
            this.useCustomSlidingRotationalSpeed =
              options.useCustomSlidingRotationalSpeed;
            this.sliding = false;
            this.chassisConnectionPointLocal =
              options.chassisConnectionPointLocal.clone();
            this.chassisConnectionPointWorld =
              options.chassisConnectionPointWorld.clone();
            this.directionLocal = options.directionLocal.clone();
            this.directionWorld = options.directionWorld.clone();
            this.axleLocal = options.axleLocal.clone();
            this.axleWorld = options.axleWorld.clone();
            this.suspensionRestLength = options.suspensionRestLength;
            this.suspensionMaxLength = options.suspensionMaxLength;
            this.radius = options.radius;
            this.suspensionStiffness = options.suspensionStiffness;
            this.dampingCompression = options.dampingCompression;
            this.dampingRelaxation = options.dampingRelaxation;
            this.frictionSlip = options.frictionSlip;
            this.forwardAcceleration = options.forwardAcceleration;
            this.sideAcceleration = options.sideAcceleration;
            this.steering = 0;
            this.rotation = 0;
            this.deltaRotation = 0;
            this.rollInfluence = options.rollInfluence;
            this.maxSuspensionForce = options.maxSuspensionForce;
            this.engineForce = 0;
            this.brake = 0;
            this.isFrontWheel = options.isFrontWheel;
            this.clippedInvContactDotSuspension = 1;
            this.suspensionRelativeVelocity = 0;
            this.suspensionForce = 0;
            this.slipInfo = 0;
            this.skidInfo = 0;
            this.suspensionLength = 0;
            this.sideImpulse = 0;
            this.forwardImpulse = 0;
            this.raycastResult = new RaycastResult();
            this.worldTransform = new Transform();
            this.isInContact = false;
          }

          updateWheel(chassis) {
            const raycastResult = this.raycastResult;

            if (this.isInContact) {
              const project = raycastResult.hitNormalWorld.dot(
                raycastResult.directionWorld
              );
              raycastResult.hitPointWorld.vsub(chassis.position, relpos);
              chassis.getVelocityAtWorldPoint(
                relpos,
                chassis_velocity_at_contactPoint
              );
              const projVel = raycastResult.hitNormalWorld.dot(
                chassis_velocity_at_contactPoint
              );

              if (project >= -0.1) {
                this.suspensionRelativeVelocity = 0.0;
                this.clippedInvContactDotSuspension = 1.0 / 0.1;
              } else {
                const inv = -1 / project;
                this.suspensionRelativeVelocity = projVel * inv;
                this.clippedInvContactDotSuspension = inv;
              }
            } else {
              // Not in contact : position wheel in a nice (rest length) position
              raycastResult.suspensionLength = this.suspensionRestLength;
              this.suspensionRelativeVelocity = 0.0;
              raycastResult.directionWorld.scale(
                -1,
                raycastResult.hitNormalWorld
              );
              this.clippedInvContactDotSuspension = 1.0;
            }
          }
        }
        const chassis_velocity_at_contactPoint = new Vec3();
        const relpos = new Vec3();

        /**
         * Vehicle helper class that casts rays from the wheel positions towards the ground and applies forces.
         */
        class RaycastVehicle {
          /** The car chassis body. */

          /** The wheels. */

          /** Will be set to true if the car is sliding. */

          /** Index of the right axis. x=0, y=1, z=2 */

          /** Index of the forward axis. x=0, y=1, z=2 */

          /** Index of the up axis. x=0, y=1, z=2 */

          /** The constraints. */

          /** Optional pre-step callback. */

          /** Number of wheels on the ground. */
          constructor(options) {
            this.chassisBody = void 0;
            this.wheelInfos = void 0;
            this.sliding = void 0;
            this.world = void 0;
            this.indexRightAxis = void 0;
            this.indexForwardAxis = void 0;
            this.indexUpAxis = void 0;
            this.constraints = void 0;
            this.preStepCallback = void 0;
            this.currentVehicleSpeedKmHour = void 0;
            this.numWheelsOnGround = void 0;
            this.chassisBody = options.chassisBody;
            this.wheelInfos = [];
            this.sliding = false;
            this.world = null;
            this.indexRightAxis =
              typeof options.indexRightAxis !== "undefined"
                ? options.indexRightAxis
                : 2;
            this.indexForwardAxis =
              typeof options.indexForwardAxis !== "undefined"
                ? options.indexForwardAxis
                : 0;
            this.indexUpAxis =
              typeof options.indexUpAxis !== "undefined"
                ? options.indexUpAxis
                : 1;
            this.constraints = [];

            this.preStepCallback = () => {};

            this.currentVehicleSpeedKmHour = 0;
            this.numWheelsOnGround = 0;
          }
          /**
           * Add a wheel. For information about the options, see `WheelInfo`.
           */

          addWheel(options = {}) {
            const info = new WheelInfo(options);
            const index = this.wheelInfos.length;
            this.wheelInfos.push(info);
            return index;
          }
          /**
           * Set the steering value of a wheel.
           */

          setSteeringValue(value, wheelIndex) {
            const wheel = this.wheelInfos[wheelIndex];
            wheel.steering = value;
          }
          /**
           * Set the wheel force to apply on one of the wheels each time step
           */

          applyEngineForce(value, wheelIndex) {
            this.wheelInfos[wheelIndex].engineForce = value;
          }
          /**
           * Set the braking force of a wheel
           */

          setBrake(brake, wheelIndex) {
            this.wheelInfos[wheelIndex].brake = brake;
          }
          /**
           * Add the vehicle including its constraints to the world.
           */

          addToWorld(world) {
            world.addBody(this.chassisBody);
            const that = this;

            this.preStepCallback = () => {
              that.updateVehicle(world.dt);
            };

            world.addEventListener("preStep", this.preStepCallback);
            this.world = world;
          }
          /**
           * Get one of the wheel axles, world-oriented.
           */

          getVehicleAxisWorld(axisIndex, result) {
            result.set(
              axisIndex === 0 ? 1 : 0,
              axisIndex === 1 ? 1 : 0,
              axisIndex === 2 ? 1 : 0
            );
            this.chassisBody.vectorToWorldFrame(result, result);
          }

          updateVehicle(timeStep) {
            const wheelInfos = this.wheelInfos;
            const numWheels = wheelInfos.length;
            const chassisBody = this.chassisBody;

            for (let i = 0; i < numWheels; i++) {
              this.updateWheelTransform(i);
            }

            this.currentVehicleSpeedKmHour =
              3.6 * chassisBody.velocity.length();
            const forwardWorld = new Vec3();
            this.getVehicleAxisWorld(this.indexForwardAxis, forwardWorld);

            if (forwardWorld.dot(chassisBody.velocity) < 0) {
              this.currentVehicleSpeedKmHour *= -1;
            } // simulate suspension

            for (let i = 0; i < numWheels; i++) {
              this.castRay(wheelInfos[i]);
            }

            this.updateSuspension(timeStep);
            const impulse = new Vec3();
            const relpos = new Vec3();

            for (let i = 0; i < numWheels; i++) {
              //apply suspension force
              const wheel = wheelInfos[i];
              let suspensionForce = wheel.suspensionForce;

              if (suspensionForce > wheel.maxSuspensionForce) {
                suspensionForce = wheel.maxSuspensionForce;
              }

              wheel.raycastResult.hitNormalWorld.scale(
                suspensionForce * timeStep,
                impulse
              );
              wheel.raycastResult.hitPointWorld.vsub(
                chassisBody.position,
                relpos
              );
              chassisBody.applyImpulse(impulse, relpos);
            }

            this.updateFriction(timeStep);
            const hitNormalWorldScaledWithProj = new Vec3();
            const fwd = new Vec3();
            const vel = new Vec3();

            for (let i = 0; i < numWheels; i++) {
              const wheel = wheelInfos[i]; //const relpos = new Vec3();
              //wheel.chassisConnectionPointWorld.vsub(chassisBody.position, relpos);

              chassisBody.getVelocityAtWorldPoint(
                wheel.chassisConnectionPointWorld,
                vel
              ); // Hack to get the rotation in the correct direction

              let m = 1;

              switch (this.indexUpAxis) {
                case 1:
                  m = -1;
                  break;
              }

              if (wheel.isInContact) {
                this.getVehicleAxisWorld(this.indexForwardAxis, fwd);
                const proj = fwd.dot(wheel.raycastResult.hitNormalWorld);
                wheel.raycastResult.hitNormalWorld.scale(
                  proj,
                  hitNormalWorldScaledWithProj
                );
                fwd.vsub(hitNormalWorldScaledWithProj, fwd);
                const proj2 = fwd.dot(vel);
                wheel.deltaRotation = (m * proj2 * timeStep) / wheel.radius;
              }

              if (
                (wheel.sliding || !wheel.isInContact) &&
                wheel.engineForce !== 0 &&
                wheel.useCustomSlidingRotationalSpeed
              ) {
                // Apply custom rotation when accelerating and sliding
                wheel.deltaRotation =
                  (wheel.engineForce > 0 ? 1 : -1) *
                  wheel.customSlidingRotationalSpeed *
                  timeStep;
              } // Lock wheels

              if (Math.abs(wheel.brake) > Math.abs(wheel.engineForce)) {
                wheel.deltaRotation = 0;
              }

              wheel.rotation += wheel.deltaRotation; // Use the old value

              wheel.deltaRotation *= 0.99; // damping of rotation when not in contact
            }
          }

          updateSuspension(deltaTime) {
            const chassisBody = this.chassisBody;
            const chassisMass = chassisBody.mass;
            const wheelInfos = this.wheelInfos;
            const numWheels = wheelInfos.length;

            for (let w_it = 0; w_it < numWheels; w_it++) {
              const wheel = wheelInfos[w_it];

              if (wheel.isInContact) {
                let force; // Spring

                const susp_length = wheel.suspensionRestLength;
                const current_length = wheel.suspensionLength;
                const length_diff = susp_length - current_length;
                force =
                  wheel.suspensionStiffness *
                  length_diff *
                  wheel.clippedInvContactDotSuspension; // Damper

                const projected_rel_vel = wheel.suspensionRelativeVelocity;
                let susp_damping;

                if (projected_rel_vel < 0) {
                  susp_damping = wheel.dampingCompression;
                } else {
                  susp_damping = wheel.dampingRelaxation;
                }

                force -= susp_damping * projected_rel_vel;
                wheel.suspensionForce = force * chassisMass;

                if (wheel.suspensionForce < 0) {
                  wheel.suspensionForce = 0;
                }
              } else {
                wheel.suspensionForce = 0;
              }
            }
          }
          /**
           * Remove the vehicle including its constraints from the world.
           */

          removeFromWorld(world) {
            this.constraints;
            world.removeBody(this.chassisBody);
            world.removeEventListener("preStep", this.preStepCallback);
            this.world = null;
          }

          castRay(wheel) {
            const rayvector = castRay_rayvector;
            const target = castRay_target;
            this.updateWheelTransformWorld(wheel);
            const chassisBody = this.chassisBody;
            let depth = -1;
            const raylen = wheel.suspensionRestLength + wheel.radius;
            wheel.directionWorld.scale(raylen, rayvector);
            const source = wheel.chassisConnectionPointWorld;
            source.vadd(rayvector, target);
            const raycastResult = wheel.raycastResult;
            raycastResult.reset(); // Turn off ray collision with the chassis temporarily

            const oldState = chassisBody.collisionResponse;
            chassisBody.collisionResponse = false; // Cast ray against world

            this.world.rayTest(source, target, raycastResult);
            chassisBody.collisionResponse = oldState;
            const object = raycastResult.body;
            wheel.raycastResult.groundObject = 0;

            if (object) {
              depth = raycastResult.distance;
              wheel.raycastResult.hitNormalWorld = raycastResult.hitNormalWorld;
              wheel.isInContact = true;
              const hitDistance = raycastResult.distance;
              wheel.suspensionLength = hitDistance - wheel.radius; // clamp on max suspension travel

              const minSuspensionLength =
                wheel.suspensionRestLength - wheel.maxSuspensionTravel;
              const maxSuspensionLength =
                wheel.suspensionRestLength + wheel.maxSuspensionTravel;

              if (wheel.suspensionLength < minSuspensionLength) {
                wheel.suspensionLength = minSuspensionLength;
              }

              if (wheel.suspensionLength > maxSuspensionLength) {
                wheel.suspensionLength = maxSuspensionLength;
                wheel.raycastResult.reset();
              }

              const denominator = wheel.raycastResult.hitNormalWorld.dot(
                wheel.directionWorld
              );
              const chassis_velocity_at_contactPoint = new Vec3();
              chassisBody.getVelocityAtWorldPoint(
                wheel.raycastResult.hitPointWorld,
                chassis_velocity_at_contactPoint
              );
              const projVel = wheel.raycastResult.hitNormalWorld.dot(
                chassis_velocity_at_contactPoint
              );

              if (denominator >= -0.1) {
                wheel.suspensionRelativeVelocity = 0;
                wheel.clippedInvContactDotSuspension = 1 / 0.1;
              } else {
                const inv = -1 / denominator;
                wheel.suspensionRelativeVelocity = projVel * inv;
                wheel.clippedInvContactDotSuspension = inv;
              }
            } else {
              //put wheel info as in rest position
              wheel.suspensionLength =
                wheel.suspensionRestLength + 0 * wheel.maxSuspensionTravel;
              wheel.suspensionRelativeVelocity = 0.0;
              wheel.directionWorld.scale(
                -1,
                wheel.raycastResult.hitNormalWorld
              );
              wheel.clippedInvContactDotSuspension = 1.0;
            }

            return depth;
          }

          updateWheelTransformWorld(wheel) {
            wheel.isInContact = false;
            const chassisBody = this.chassisBody;
            chassisBody.pointToWorldFrame(
              wheel.chassisConnectionPointLocal,
              wheel.chassisConnectionPointWorld
            );
            chassisBody.vectorToWorldFrame(
              wheel.directionLocal,
              wheel.directionWorld
            );
            chassisBody.vectorToWorldFrame(wheel.axleLocal, wheel.axleWorld);
          }
          /**
           * Update one of the wheel transform.
           * Note when rendering wheels: during each step, wheel transforms are updated BEFORE the chassis; ie. their position becomes invalid after the step. Thus when you render wheels, you must update wheel transforms before rendering them. See raycastVehicle demo for an example.
           * @param wheelIndex The wheel index to update.
           */

          updateWheelTransform(wheelIndex) {
            const up = tmpVec4;
            const right = tmpVec5;
            const fwd = tmpVec6;
            const wheel = this.wheelInfos[wheelIndex];
            this.updateWheelTransformWorld(wheel);
            wheel.directionLocal.scale(-1, up);
            right.copy(wheel.axleLocal);
            up.cross(right, fwd);
            fwd.normalize();
            right.normalize(); // Rotate around steering over the wheelAxle

            const steering = wheel.steering;
            const steeringOrn = new Quaternion();
            steeringOrn.setFromAxisAngle(up, steering);
       