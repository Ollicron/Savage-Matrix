#[
The purpose of this module is to implement 3D and 4D matrices, vectors are represented in row order, meaning one row can represent one vector.
]#

#[ Imports ]#
import vector_types

type
  Matrix3D* = object
    mat*: array[3, array[3, float]]
  Matrix4D* = object
    mat*: array[4, array[4, float]]
  Transform* = object
    mat*: array[4, array[4, float]]

func matrix3D*(A: array[3, float|int], B: array[3, float|int], C: array[3,
    float|int]): Matrix3D {.inline.} =
  result.mat[0][0] = float(A[0]); result.mat[0][1] = float(A[1]); result.mat[0][
      2] = float(A[2]);
  result.mat[1][0] = float(B[0]); result.mat[1][1] = float(B[1]); result.mat[1][
      2] = float(B[2]);
  result.mat[2][0] = float(C[0]); result.mat[2][1] = float(C[1]); result.mat[2][
      2] = float(C[2]);

func matrix3D*(A: Vector3D, B: Vector3D, C: Vector3D): Matrix3D {.inline.} =
  result.mat[0][0] = A.x; result.mat[0][1] = A.y; result.mat[0][2] = A.z;
  result.mat[1][0] = B.x; result.mat[1][1] = B.y; result.mat[1][2] = B.z;
  result.mat[2][0] = C.x; result.mat[2][1] = C.y; result.mat[2][2] = C.z;

# Constructors for Matrix4D
func matrix4D*(A, B, C, D: array[4, float|int]): Matrix4D {.inline.} =
  # Row 0
  result.mat[0][0] = float(A[0]); result.mat[0][1] = float(A[1])
  result.mat[0][2] = float(A[2]); result.mat[0][3] = float(A[3])
  # Row 1
  result.mat[1][0] = float(B[0]); result.mat[1][1] = float(B[1])
  result.mat[1][2] = float(B[2]); result.mat[1][3] = float(B[3])
  # Row 2
  result.mat[2][0] = float(C[0]); result.mat[2][1] = float(C[1])
  result.mat[2][2] = float(C[2]); result.mat[2][3] = float(C[3])
  # Row 3
  result.mat[3][0] = float(D[0]); result.mat[3][1] = float(D[1])
  result.mat[3][2] = float(D[2]); result.mat[3][3] = float(D[3])

func matrix4D*(A, B, C, D: Vector4D): Matrix4D {.inline.} =
  # Row 0
  result.mat[0][0] = A.x; result.mat[0][1] = A.y
  result.mat[0][2] = A.z; result.mat[0][3] = A.w
  # Row 1
  result.mat[1][0] = B.x; result.mat[1][1] = B.y
  result.mat[1][2] = B.z; result.mat[1][3] = B.w
  # Row 2
  result.mat[2][0] = C.x; result.mat[2][1] = C.y
  result.mat[2][2] = C.z; result.mat[2][3] = C.w
  # Row 3
  result.mat[3][0] = D.x; result.mat[3][1] = D.y
  result.mat[3][2] = D.z; result.mat[3][3] = D.w


# 3x3 direct constructor
func matrix3D*(n00, n01, n02,
               n10, n11, n12,
               n20, n21, n22: float): Matrix3D =
  result.mat[0][0] = n00
  result.mat[0][1] = n01
  result.mat[0][2] = n02
  result.mat[1][0] = n10
  result.mat[1][1] = n11
  result.mat[1][2] = n12
  result.mat[2][0] = n20
  result.mat[2][1] = n21
  result.mat[2][2] = n22

# 4x4 direct constructor
func matrix4D*(n00, n01, n02, n03,
               n10, n11, n12, n13,
               n20, n21, n22, n23,
               n30, n31, n32, n33: float): Matrix4D =
  result.mat[0][0] = n00
  result.mat[0][1] = n01
  result.mat[0][2] = n02
  result.mat[0][3] = n03
  result.mat[1][0] = n10
  result.mat[1][1] = n11
  result.mat[1][2] = n12
  result.mat[1][3] = n13
  result.mat[2][0] = n20
  result.mat[2][1] = n21
  result.mat[2][2] = n22
  result.mat[2][3] = n23
  result.mat[3][0] = n30
  result.mat[3][1] = n31
  result.mat[3][2] = n32
  result.mat[3][3] = n33


# Constructor for a transformation matrix using a translation vector and 3D matrix. The convention is row order so tvec is at the bottom of the matrix
func transform*(rotation:Matrix3D,tvec:Vector3D):Transform {.inline.}=
  result.mat[0][0] = rotation.mat[0][0]
  result.mat[0][1] = rotation.mat[0][1]
  result.mat[0][2] = rotation.mat[0][2]
  result.mat[1][0] = rotation.mat[1][0]
  result.mat[1][1] = rotation.mat[1][1]
  result.mat[1][2] = rotation.mat[1][2]
  result.mat[2][0] = rotation.mat[2][0]
  result.mat[2][1] = rotation.mat[2][1]
  result.mat[2][2] = rotation.mat[2][2]
  result.mat[3][0] = tvec[0]
  result.mat[3][1] = tvec[1]
  result.mat[3][2] = tvec[2]
  result.mat[3][3] = 1.0
  result.mat[0][3] = 0.0
  result.mat[1][3] = 0.0
  result.mat[2][3] = 0.0

# operators for accessing and mutating elements for transforms
func`[]`*(transform: Transform,A:int, B: int):float {.inline.}=
  result = transform.mat[A][B]

func`[]=`*(transform: var Transform,A:int,B: int, val:float) {.inline.} =
  transform.mat[A][B]=val


# Access operator
{.experimental: "callOperator".}
func `()`*(matrix: Matrix4D, A: int, B: int): float {.inline.} =
  result = matrix.mat[A][B]

{.experimental: "callOperator".}
func `()`*(matrix: Matrix3D, A: int, B: int): float {.inline.} =
  result = matrix.mat[A][B]

when isMainModule:
  let matr = Matrix3D(mat: [[1, 2, 3], [4, 5, 6], [7, 8, 9]])
  echo matr.mat[1]
  echo matr(1, 1)
