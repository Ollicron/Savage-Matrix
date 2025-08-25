#[
The purpose of this module is to implement the mathematical operations for the matrix_types module. It is meant to be used for fast calculations. Skewing operations are a scam. Do not use them. All rotations are performed using the left hand rule.
]#



#[ Imports ]#
import
  matrix_types,
  vector_types,
  vector_math,
  std/[math]



#[ Templates ]#




# Template to make basic matrix operations
template makeMatrixOps(op) =
  func op(A: Matrix3D, B: Matrix3D): Matrix3D {.inline.} =
    result.mat[0][0] = op(A(0, 0), B(0, 0))
    result.mat[0][1] = op(A(0, 1), B(0, 1))
    result.mat[0][2] = op(A(0, 2), B(0, 2))
    result.mat[1][0] = op(A(1, 0), B(1, 0))
    result.mat[1][1] = op(A(1, 1), B(1, 1))
    result.mat[1][2] = op(A(1, 2), B(1, 2))
    result.mat[2][0] = op(A(2, 0), B(2, 0))
    result.mat[2][1] = op(A(2, 1), B(2, 1))
    result.mat[2][2] = op(A(2, 2), B(2, 2))

makeMatrixOps(`+`)
makeMatrixOps(`-`)




#[ Procedures/Functions ]#




# Function to scale a matrix
func scale(A: Matrix3D, val: float|int): Matrix3D{.inline.} =
  result.mat[0][0] = A(0, 0) * float(val)
  result.mat[0][1] = A(0, 1) * float(val)
  result.mat[0][2] = A(0, 2) * float(val)
  result.mat[1][0] = A(1, 0) * float(val)
  result.mat[1][1] = A(1, 1) * float(val)
  result.mat[1][2] = A(1, 2) * float(val)
  result.mat[2][0] = A(2, 0) * float(val)
  result.mat[2][1] = A(2, 1) * float(val)
  result.mat[2][2] = A(2, 2) * float(val)

func scale(A: Matrix4D, val: float|int): Matrix4D{.inline.} =
  result.mat[0][0] = A(0, 0) * float(val)
  result.mat[0][1] = A(0, 1) * float(val)
  result.mat[0][2] = A(0, 2) * float(val)
  result.mat[0][3] = A(0, 3) * float(val)
  result.mat[1][0] = A(1, 0) * float(val)
  result.mat[1][1] = A(1, 1) * float(val)
  result.mat[1][2] = A(1, 2) * float(val)
  result.mat[1][3] = A(1, 3) * float(val)
  result.mat[2][0] = A(2, 0) * float(val)
  result.mat[2][1] = A(2, 1) * float(val)
  result.mat[2][2] = A(2, 2) * float(val)
  result.mat[2][3] = A(2, 3) * float(val)
  result.mat[3][0] = A(3, 0) * float(val)
  result.mat[3][1] = A(3, 1) * float(val)
  result.mat[3][2] = A(3, 2) * float(val)
  result.mat[3][3] = A(3, 3) * float(val)

# Function to multiply matrices
func `*`(A: Matrix3D, B: Matrix3D): Matrix3D {.inline.} =
  # Row 0
  result.mat[0][0] = A(0, 0)*B(0, 0) + A(0, 1)*B(1, 0) + A(0, 2)*B(2, 0)
  result.mat[0][1] = A(0, 0)*B(0, 1) + A(0, 1)*B(1, 1) + A(0, 2)*B(2, 1)
  result.mat[0][2] = A(0, 0)*B(0, 2) + A(0, 1)*B(1, 2) + A(0, 2)*B(2, 2)

  # Row 1
  result.mat[1][0] = A(1, 0)*B(0, 0) + A(1, 1)*B(1, 0) + A(1, 2)*B(2, 0)
  result.mat[1][1] = A(1, 0)*B(0, 1) + A(1, 1)*B(1, 1) + A(1, 2)*B(2, 1)
  result.mat[1][2] = A(1, 0)*B(0, 2) + A(1, 1)*B(1, 2) + A(1, 2)*B(2, 2)

  # Row 2
  result.mat[2][0] = A(2, 0)*B(0, 0) + A(2, 1)*B(1, 0) + A(2, 2)*B(2, 0)
  result.mat[2][1] = A(2, 0)*B(0, 1) + A(2, 1)*B(1, 1) + A(2, 2)*B(2, 1)
  result.mat[2][2] = A(2, 0)*B(0, 2) + A(2, 1)*B(1, 2) + A(2, 2)*B(2, 2)

func `*`(A, B: Matrix4D): Matrix4D {.inline.} =
  # Row 0
  result.mat[0][0] = A.mat[0][0]*B.mat[0][0] + A.mat[0][1]*B.mat[1][0] + A.mat[
      0][2]*B.mat[2][0] + A.mat[0][3]*B.mat[3][0]
  result.mat[0][1] = A.mat[0][0]*B.mat[0][1] + A.mat[0][1]*B.mat[1][1] + A.mat[
      0][2]*B.mat[2][1] + A.mat[0][3]*B.mat[3][1]
  result.mat[0][2] = A.mat[0][0]*B.mat[0][2] + A.mat[0][1]*B.mat[1][2] + A.mat[
      0][2]*B.mat[2][2] + A.mat[0][3]*B.mat[3][2]
  result.mat[0][3] = A.mat[0][0]*B.mat[0][3] + A.mat[0][1]*B.mat[1][3] + A.mat[
      0][2]*B.mat[2][3] + A.mat[0][3]*B.mat[3][3]

  # Row 1
  result.mat[1][0] = A.mat[1][0]*B.mat[0][0] + A.mat[1][1]*B.mat[1][0] + A.mat[
      1][2]*B.mat[2][0] + A.mat[1][3]*B.mat[3][0]
  result.mat[1][1] = A.mat[1][0]*B.mat[0][1] + A.mat[1][1]*B.mat[1][1] + A.mat[
      1][2]*B.mat[2][1] + A.mat[1][3]*B.mat[3][1]
  result.mat[1][2] = A.mat[1][0]*B.mat[0][2] + A.mat[1][1]*B.mat[1][2] + A.mat[
      1][2]*B.mat[2][2] + A.mat[1][3]*B.mat[3][2]
  result.mat[1][3] = A.mat[1][0]*B.mat[0][3] + A.mat[1][1]*B.mat[1][3] + A.mat[
      1][2]*B.mat[2][3] + A.mat[1][3]*B.mat[3][3]

  # Row 2
  result.mat[2][0] = A.mat[2][0]*B.mat[0][0] + A.mat[2][1]*B.mat[1][0] + A.mat[
      2][2]*B.mat[2][0] + A.mat[2][3]*B.mat[3][0]
  result.mat[2][1] = A.mat[2][0]*B.mat[0][1] + A.mat[2][1]*B.mat[1][1] + A.mat[
      2][2]*B.mat[2][1] + A.mat[2][3]*B.mat[3][1]
  result.mat[2][2] = A.mat[2][0]*B.mat[0][2] + A.mat[2][1]*B.mat[1][2] + A.mat[
      2][2]*B.mat[2][2] + A.mat[2][3]*B.mat[3][2]
  result.mat[2][3] = A.mat[2][0]*B.mat[0][3] + A.mat[2][1]*B.mat[1][3] + A.mat[
      2][2]*B.mat[2][3] + A.mat[2][3]*B.mat[3][3]

  # Row 3
  result.mat[3][0] = A.mat[3][0]*B.mat[0][0] + A.mat[3][1]*B.mat[1][0] + A.mat[
      3][2]*B.mat[2][0] + A.mat[3][3]*B.mat[3][0]
  result.mat[3][1] = A.mat[3][0]*B.mat[0][1] + A.mat[3][1]*B.mat[1][1] + A.mat[
      3][2]*B.mat[2][1] + A.mat[3][3]*B.mat[3][1]
  result.mat[3][2] = A.mat[3][0]*B.mat[0][2] + A.mat[3][1]*B.mat[1][2] + A.mat[
      3][2]*B.mat[2][2] + A.mat[3][3]*B.mat[3][2]
  result.mat[3][3] = A.mat[3][0]*B.mat[0][3] + A.mat[3][1]*B.mat[1][3] + A.mat[
      3][2]*B.mat[2][3] + A.mat[3][3]*B.mat[3][3]

# Function to multiply Matrices to vectors.
func `*`(A: Matrix3D, B: Vector3D): Vector3D {.inline.} =
  result.x = A(0, 0)*B.x + A(0, 1)*B.y + A(0, 2)*B.z
  result.y = A(1, 0)*B.x + A(1, 1)*B.y + A(1, 2)*B.z
  result.z = A(2, 0)*B.x + A(2, 1)*B.y + A(2, 2)*B.z

func `*`(A: Matrix4D, B: Vector4D): Vector4D {.inline.} =
  result.x = A(0, 0)*B.x + A(0, 1)*B.y + A(0, 2)*B.z + A(0, 3)*B.w
  result.y = A(1, 0)*B.x + A(1, 1)*B.y + A(1, 2)*B.z + A(1, 3)*B.w
  result.z = A(2, 0)*B.x + A(2, 1)*B.y + A(2, 2)*B.z + A(2, 3)*B.w
  result.w = A(3, 0)*B.x + A(3, 1)*B.y + A(3, 2)*B.z + A(3, 3)*B.w

# Function for the determinant of a 3x3 matrix
func det*(M: Matrix3D): float {.inline.} =
  let a = M.mat
  a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1]) -
  a[0][1] * (a[1][0]*a[2][2] - a[1][2]*a[2][0]) +
  a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0])

# Function to compute determinant of a 4x4 matrix
func det*(M: Matrix4D): float {.inline.} =
  let m = M.mat
  # Assign elements for readability
  let m00 = m[0][0].float
  let m01 = m[0][1].float
  let m02 = m[0][2].float
  let m03 = m[0][3].float

  let m10 = m[1][0].float
  let m11 = m[1][1].float
  let m12 = m[1][2].float
  let m13 = m[1][3].float

  let m20 = m[2][0].float
  let m21 = m[2][1].float
  let m22 = m[2][2].float
  let m23 = m[2][3].float

  let m30 = m[3][0].float
  let m31 = m[3][1].float
  let m32 = m[3][2].float
  let m33 = m[3][3].float

  # Compute determinant via expansion by first row
  result =
    m00 * (m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 -
        m22*m31)) -
    m01 * (m10*(m22*m33 - m23*m32) - m12*(m20*m33 - m23*m30) + m13*(m20*m32 -
        m22*m30)) +
    m02 * (m10*(m21*m33 - m23*m31) - m11*(m20*m33 - m23*m30) + m13*(m20*m31 -
        m21*m30)) -
    m03 * (m10*(m21*m32 - m22*m31) - m11*(m20*m32 - m22*m30) + m12*(m20*m31 - m21*m30))

# Function to transpose a matrix
func transpose*(M: Matrix3D): Matrix3D {.inline.}=
  matrix3D(
    vector(M.mat[0][0], M.mat[1][0], M.mat[2][0]),
    vector(M.mat[0][1], M.mat[1][1], M.mat[2][1]),
    vector(M.mat[0][2], M.mat[1][2], M.mat[2][2])
  )

func transpose*(M: Matrix4D): Matrix4D {.inline.}=
  matrix4d(
    vector(M.mat[0][0], M.mat[1][0], M.mat[2][0], M.mat[3][0]),
    vector(M.mat[0][1], M.mat[1][1], M.mat[2][1], M.mat[3][1]),
    vector(M.mat[0][2], M.mat[1][2], M.mat[2][2], M.mat[3][2]),
    vector(M.mat[0][3], M.mat[1][3], M.mat[2][3], M.mat[3][3])
  )

# Function to get the inverse of a matrix (FAST direct access only)
func inv*(M: Matrix3D): Matrix3D {.inline.} =
  let m = M.mat

  let m00 = m[0][0].float
  let m01 = m[0][1].float
  let m02 = m[0][2].float
  let m10 = m[1][0].float
  let m11 = m[1][1].float
  let m12 = m[1][2].float
  let m20 = m[2][0].float
  let m21 = m[2][1].float
  let m22 = m[2][2].float

  let det = m00*(m11*m22 - m12*m21) -
            m01*(m10*m22 - m12*m20) +
            m02*(m10*m21 - m11*m20)

  let invDet = 1.0 / det

  result = matrix3d(
    vector((m11*m22 - m12*m21) * invDet,
           -(m01*m22 - m02*m21) * invDet,
           (m01*m12 - m02*m11) * invDet),
    vector(-(m10*m22 - m12*m20) * invDet,
           (m00*m22 - m02*m20) * invDet,
           -(m00*m12 - m02*m10) * invDet),
    vector((m10*m21 - m11*m20) * invDet,
           -(m00*m21 - m01*m20) * invDet,
           (m00*m11 - m01*m10) * invDet)
  )

# Inverse of a 4x4 matrix using the determinant function det4
func inv*(M: Matrix4D): Matrix4D =
  let m = M.mat
  let det = det(M)
  if abs(det) < 1e-8:
    raise newException(ValueError, "Matrix is singular and cannot be inverted")
  let invDet = 1.0 / det

  # Precompute all the 2x2 minors needed for cofactors (first two rows + last two)
  let c00 = m[2][2]*m[3][3] - m[2][3]*m[3][2]
  let c01 = m[2][1]*m[3][3] - m[2][3]*m[3][1]
  let c02 = m[2][1]*m[3][2] - m[2][2]*m[3][1]
  let c03 = m[2][0]*m[3][3] - m[2][3]*m[3][0]
  let c04 = m[2][0]*m[3][2] - m[2][2]*m[3][0]
  let c05 = m[2][0]*m[3][1] - m[2][1]*m[3][0]

  # Build the inverse matrix directly using the 4x4 constructor
  result = matrix4D(
    (m[1][1]*c00 - m[1][2]*c01 + m[1][3]*c02) * invDet,
    (-m[0][1]*c00 + m[0][2]*c01 - m[0][3]*c02) * invDet,
    (m[0][1]*(m[1][2]*m[3][3] - m[1][3]*m[3][2]) - m[0][2]*(m[1][1]*m[3][3] - m[
        1][3]*m[3][1]) + m[0][3]*(m[1][1]*m[3][2] - m[1][2]*m[3][1])) * invDet,
    (-m[0][1]*(m[1][2]*m[2][3] - m[1][3]*m[2][2]) + m[0][2]*(m[1][1]*m[2][3] -
        m[1][3]*m[2][1]) - m[0][3]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])) *
        invDet,

    (-m[1][0]*c00 + m[1][2]*c03 - m[1][3]*c04) * invDet,
    (m[0][0]*c00 - m[0][2]*c03 + m[0][3]*c04) * invDet,
    (-m[0][0]*(m[1][2]*m[3][3] - m[1][3]*m[3][2]) + m[0][2]*(m[1][0]*m[3][3] -
        m[1][3]*m[3][0]) - m[0][3]*(m[1][0]*m[3][2] - m[1][2]*m[3][0])) *
        invDet,
    (m[0][0]*(m[1][2]*m[2][3] - m[1][3]*m[2][2]) - m[0][2]*(m[1][0]*m[2][3] - m[
        1][3]*m[2][0]) + m[0][3]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])) * invDet,

    (m[1][0]*c01 - m[1][1]*c03 + m[1][3]*c05) * invDet,
    (-m[0][0]*c01 + m[0][1]*c03 - m[0][3]*c05) * invDet,
    (m[0][0]*(m[1][1]*m[3][3] - m[1][3]*m[3][1]) - m[0][1]*(m[1][0]*m[3][3] - m[
        1][3]*m[3][0]) + m[0][3]*(m[1][0]*m[3][1] - m[1][1]*m[3][0])) * invDet,
    (-m[0][0]*(m[1][1]*m[2][3] - m[1][3]*m[2][1]) + m[0][1]*(m[1][0]*m[2][3] -
        m[1][3]*m[2][0]) - m[0][3]*(m[1][0]*m[2][1] - m[1][1]*m[2][0])) *
        invDet,

    (-m[1][0]*c02 + m[1][1]*c04 - m[1][2]*c05) * invDet,
    (m[0][0]*c02 - m[0][1]*c04 + m[0][2]*c05) * invDet,
    (-m[0][0]*(m[1][1]*m[3][2] - m[1][2]*m[3][1]) + m[0][1]*(m[1][0]*m[3][2] -
        m[1][2]*m[3][0]) - m[0][2]*(m[1][0]*m[3][1] - m[1][1]*m[3][0])) *
        invDet,
    (m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2] - m[
        1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0])) * invDet
  )


# We will only ever need the 3D versions of rotation matrices
func makeRotX*(A: float): Matrix3D {.inline.} =
  let c = cos(A)
  let s = sin(A)
  return matrix3D(
      1.0, 0.0, 0.0,
      0.0, c, s,   # flip the sign here for left-hand
      0.0, -s, c
  )

func makeRotY*(A: float): Matrix3D {.inline.} =
  let c = cos(A)
  let s = sin(A)
  return matrix3D(
      c, 0.0, -s,  # flip sign
      0.0, 1.0, 0.0,
      s, 0.0, c
  )

func makeRotZ*(A: float): Matrix3D {.inline.} =
  let c = cos(A)
  let s = sin(A)
  return matrix3D(
      c, s, 0.0,   # flip sign
      -s, c, 0.0,
      0.0, 0.0, 1.0
  )

func rotOnAxis*(unitAxis: Vector3D, theta: float): Matrix3D {.inline.}=
  let n = normalize(unitAxis)
  let Nx = n.x
  let Ny = n.y
  let Nz = n.z
  let c = cos(theta)
  let s = sin(theta)
  let d = 1 - c

  result = matrix3D(
    [c + d*Nx*Nx, d*Nx*Ny + s*Nz, d*Nx*Nz - s*Ny],
    [d*Nx*Ny - s*Nz, c + d*Ny*Ny, d*Ny*Nz + s*Nx],
    [d*Nx*Nz + s*Ny, d*Ny*Nz - s*Nx, c + d*Nz*Nz]
  )


# Function to scale a basis along a specific Axis (unit vector)
func scaleOnAxis*(k: float, unitAxis: Vector3D): Matrix4D {.inline.}=
  let n = normalize(unitAxis)
  let Nx = n.x
  let Ny = n.y
  let Nz = n.z

  # Non-uniform scaling along an arbitrary axis
  result = matrix4D(
    [1 + (k-1)*Nx*Nx,     (k-1)*Nx*Ny,     (k-1)*Nx*Nz, 0.0],
    [(k-1)*Nx*Ny, 1 + (k-1)*Ny*Ny,     (k-1)*Ny*Nz, 0.0],
    [(k-1)*Nx*Nz,     (k-1)*Ny*Nz, 1 + (k-1)*Nz*Nz, 0.0],
    [0.0, 0.0, 0.0, 1.0]
  )

## Function to get matrix to  perform orthographic projection onto a plane perpendicular to unitAxis (which allows us to view object from sides)
func orthoProjectAlongAxis*(unitAxis:Vector3D): Matrix4D {.inline.} =

  let n = normalize(unitAxis)
  let Nx = n.x
  let Ny = n.y
  let Nz = n.z

  result = matrix4D(
    [1 - Nx*Nx,    -Nx*Ny,      -Nx*Nz,    0.0],
    [-Nx*Ny,       1 - Ny*Ny,   -Ny*Nz,    0.0],
    [-Nx*Nz,       -Ny*Nz,      1 - Nz*Nz, 0.0],
    [0.0,          0.0,         0.0,       1.0]
  )

# Function to get the reflection matrix from unit axis normal to the plane you want to reflect
func reflect*(axis: Vector3D): Matrix3D {.inline.}=
  let n = normalize(axis)
  let Nx = n.x
  let Ny = n.y
  let Nz = n.z

  result = matrix3D(
    [1 - 2*Nx*Nx, -2*Nx*Ny, -2*Nx*Nz],
    [-2*Nx*Ny, 1 - 2*Ny*Ny, -2*Ny*Nz],
    [-2*Nx*Nz, -2*Ny*Nz, 1 - 2*Nz*Nz]
  )



### Transformation Matrices Functions###



# Function to multiply transformation matrices
func `*`*(a, b: Transform): Transform {.inline.}=
  # row 0
  result[0,0] = a[0,0]*b[0,0] + a[0,1]*b[1,0] + a[0,2]*b[2,0] + a[0,3]*b[3,0]
  result[0,1] = a[0,0]*b[0,1] + a[0,1]*b[1,1] + a[0,2]*b[2,1] + a[0,3]*b[3,1]
  result[0,2] = a[0,0]*b[0,2] + a[0,1]*b[1,2] + a[0,2]*b[2,2] + a[0,3]*b[3,2]
  result[0,3] = a[0,0]*b[0,3] + a[0,1]*b[1,3] + a[0,2]*b[2,3] + a[0,3]*b[3,3]

  # row 1
  result[1,0] = a[1,0]*b[0,0] + a[1,1]*b[1,0] + a[1,2]*b[2,0] + a[1,3]*b[3,0]
  result[1,1] = a[1,0]*b[0,1] + a[1,1]*b[1,1] + a[1,2]*b[2,1] + a[1,3]*b[3,1]
  result[1,2] = a[1,0]*b[0,2] + a[1,1]*b[1,2] + a[1,2]*b[2,2] + a[1,3]*b[3,2]
  result[1,3] = a[1,0]*b[0,3] + a[1,1]*b[1,3] + a[1,2]*b[2,3] + a[1,3]*b[3,3]

  # row 2
  result[2,0] = a[2,0]*b[0,0] + a[2,1]*b[1,0] + a[2,2]*b[2,0] + a[2,3]*b[3,0]
  result[2,1] = a[2,0]*b[0,1] + a[2,1]*b[1,1] + a[2,2]*b[2,1] + a[2,3]*b[3,1]
  result[2,2] = a[2,0]*b[0,2] + a[2,1]*b[1,2] + a[2,2]*b[2,2] + a[2,3]*b[3,2]
  result[2,3] = a[2,0]*b[0,3] + a[2,1]*b[1,3] + a[2,2]*b[2,3] + a[2,3]*b[3,3]

  # row 3
  result[3,0] = a[3,0]*b[0,0] + a[3,1]*b[1,0] + a[3,2]*b[2,0] + a[3,3]*b[3,0]
  result[3,1] = a[3,0]*b[0,1] + a[3,1]*b[1,1] + a[3,2]*b[2,1] + a[3,3]*b[3,1]
  result[3,2] = a[3,0]*b[0,2] + a[3,1]*b[1,2] + a[3,2]*b[2,2] + a[3,3]*b[3,2]
  result[3,3] = a[3,0]*b[0,3] + a[3,1]*b[1,3] + a[3,2]*b[2,3] + a[3,3]*b[3,3]


func inverse(t: Transform): Transform{.inline.} =
  let
    a = Vector3D(x: t[0,0], y: t[0,1], z: t[0,2])
    b = Vector3D(x: t[1,0], y: t[1,1], z: t[1,2])
    c = Vector3D(x: t[2,0], y: t[2,1], z: t[2,2])
    d = Vector3D(x: t[0,3], y: t[1,3], z: t[2,3])  # translation in last column

  var s = cross(a, b)
  var tvec = cross(c, d)
  let invDet = 1.0 / dot(s, c)

  s = Vector3D(x: s.x*invDet, y: s.y*invDet, z: s.z*invDet)
  tvec = Vector3D(x: tvec.x*invDet, y: tvec.y*invDet, z: tvec.z*invDet)
  let v = Vector3D(x: c.x*invDet, y: c.y*invDet, z: c.z*invDet)

  let r0 = cross(b, v)
  let r1 = cross(v, a)

  result[0,0] = r0.x
  result[0,1] = r0.y
  result[0,2] = r0.z
  result[0,3] = -dot(b, tvec)

  result[1,0] = r1.x
  result[1,1] = r1.y
  result[1,2] = r1.z
  result[1,3] = dot(a, tvec)

  result[2,0] = s.x
  result[2,1] = s.y
  result[2,2] = s.z
  result[2,3] = -dot(d, s)

  # bottom row = [0,0,0,1] for affine
  result[3,0] = 0.0
  result[3,1] = 0.0
  result[3,2] = 0.0
  result[3,3] = 1.0

# Function to move a point using a transform
func `*`(t:Transform,p:Point3D):Point3D{.inline.}=
  result = point3D(
    t[0,0]*p.x + t[0,1]*p.y + t[0,2]*p.z + t[0,3],
    t[1,0]*p.x + t[1,1]*p.y + t[1,2]*p.z + t[1,3],
    t[2,0]*p.x + t[2,1]*p.y + t[2,2]*p.z + t[2,3]
  )

# Function to move a vector using a transform
func `*`(t:Transform,v:Vector3D):Vector3D{.inline.}=
  result = point3D(
    t[0,0]*v.x + t[0,1]*v.y + t[0,2]*v.z + t[0,3],
    t[1,0]*v.x + t[1,1]*v.y + t[1,2]*v.z + t[1,3],
    t[2,0]*v.x + t[2,1]*v.y + t[2,2]*v.z + t[2,3]
  )

# Function to get a rotation matrix from a normalized quaternion
proc getRotationMatrixLH(quat: Quat): Matrix3D =
  let q = normalizeQuat(quat)
  
  let
    x2 = q.x * q.x
    y2 = q.y * q.y
    z2 = q.z * q.z
    xy = q.x * q.y
    xz = q.x * q.z
    yz = q.y * q.z
    wx = q.w * q.x
    wy = q.w * q.y
    wz = q.w * q.z

  # Flip the signs that correspond to the left-hand rule
  result = [
    [1.0 - 2.0*(y2 + z2),  2.0*(xy + wz),      2.0*(xz - wy)],
    [2.0*(xy - wz),        1.0 - 2.0*(x2 + z2), 2.0*(yz + wx)],
    [2.0*(xz + wy),        2.0*(yz - wx),      1.0 - 2.0*(x2 + y2)]
  ]



when isMainModule:

  var H = Transform(mat:[
    [1.0, 7.0, 9.0, 0.0],  # row 0
    [2.0, 1.0, 0.0, 0.0],  # row 1
    [3.0, 4.0, 1.0, 0.0],  # row 2
    [1.0, 9.0, 3.0, 1.0]   # row 3
  ])

  let Hinv = inverse(H)
  let identity = H * Hinv

  echo "Transform * Inverse ="
  echo(identity)