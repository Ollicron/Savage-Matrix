#[
The purpose of this module is to implement behaviors and operators used for the operation of vectors. To give better performance for the library operations using loops are avoided to make direct calculations via CPU threads. Templates introduce very minor overhead only in compilation, and only in the total size of the resulting binary.
]#

#[ Imports ]#


import
  vector_types,
  std/[math],
  os, times

# Create the first set of vector ops
template makeVectorOps(op: untyped) =
  func op*(A, B: Vector2D): Vector2D {.inline.} =
    result = Vector2D(x: op(A[0], B[0]), y: op(A[1], B[1]))

  func op*(A, B: Vector3D): Vector3D {.inline.} =
    result = Vector3D(x: op(A[0], B[0]), y: op(A[1], B[1]), z: op(A[2], B[2]))

  func op*(A, B: Vector4D): Vector4D {.inline.} =
    result = Vector4D(x: op(A[0], B[0]), y: op(A[1], B[1]), z: op(A[2], B[2]),
        w: op(A[3], B[3]))



makeVectorOps(`+`)
makeVectorOps(`-`)
makeVectorOps(`*`)
makeVectorOps(`/`)
makeVectorOps(`min`)
makeVectorOps(`max`)
makeVectorOps(`mod`)

#[ Procedures/Functions ]#



## Functions to scale vectors
func scale*(A: Vector2D, val: int|float): Vector2D {.inline.} =
  result.x = A.x * float(val)
  result.y = A.y * float(val)

func scale*(A: Vector3D, val: int|float): Vector3D {.inline.} =
  result.x = A.x * float(val)
  result.y = A.y * float(val)
  result.z = A.z * float(val)

func scale*(A: Vector4D, val: int|float): Vector4D {.inline.} =
  result.x = A.x * float(val)
  result.y = A.y * float(val)
  result.z = A.z * float(val)
  result.w = A.w * float(val)


# Functions to calculate magnitude
func getMag(A: Vector2D): float {.inline.} =
  sqrt(A.x*A.x+A.y*A.y)

func getMag(A: Vector3D): float {.inline.} =
  sqrt(A.x*A.x+A.y*A.y+A.z*A.z)

func getMag(A: Vector4D): float {.inline.} =
  sqrt(A.x*A.x+A.y*A.y+A.z*A.z+A.w*A.w)


## Functions to Normalize vectors
func normalize*(A: Vector2D): Vector2D {.inline discardable.} =
  let mag = getMag(A)
  result.x = A.x / mag
  result.y = A.y / mag

## Functions to Normalize vectors
func normalize*(A: Vector3D): Vector3D {.inline discardable.} =
  let mag = getMag(A)
  result.x = A.x / mag
  result.y = A.y / mag
  result.z = A.z / mag

func normalize*(A: Vector4D): Vector4D {.inline discardable.} =
  let mag = getMag(A)
  result.x = A.x / mag
  result.y = A.y / mag
  result.z = A.z / mag
  result.w = A.w / mag


# Functions to calculate dot product between two vectors
func dot*(A, B: Vector2D): float {.inline.} =
  A.x*B.x+A.y*B.y

func dot*(A, B: Vector3D): float {.inline.} =
  A.x*B.x+A.y*B.y+A.z*B.z

func dot*(A, B: Vector4D): float {.inline.} =
  A.x*B.x+A.y*B.y+A.z*B.z+A.w*B.w

# Function to calculate left handed cross product between two vectors (only applicable to 3D vectors in graphics programming)
func cross*(A, B: Vector3D): Vector3D {.inline.} =
  result.x = A.y*B.z - A.z*B.y
  result.y = -(A.x*B.z - A.z*B.x)  # flip the sign here
  result.z = A.x*B.y - A.y*B.x

# Functions to add and remove from 3D points using vectors and other points
func `+`*(A, B: Point3D): Point3D {.inline.} =
    result = Point3D(x: `+`(A.x, B.x), y: `+`(A.y, B.y), z: `+`(A.z, B.z))

func `+`*(A: Point3D, B: Vector3D): Point3D {.inline.} =
    result = Point3D(x: `+`(A.x, B.x), y: `+`(A.y, B.y), z: `+`(A.z, B.z))

func `-`*(A, B: Point3D): Point3D {.inline.} =
    result = Point3D(x: `-`(A.x, B.x), y: `-`(A.y, B.y), z: `-`(A.z, B.z))

func `-`*(A: Point3D, B: Vector3D): Point3D {.inline.} =
    result = Point3D(x: `-`(A.x, B.x), y: `-`(A.y, B.y), z: `-`(A.z, B.z))

# Function to multiply quaternions together; the purpose is to perform rotation with q1 first then q2
func `*`*(q1, q2: Quat): Quat {.inline.} =
  result.x =  q1.w*q2.x + q1.x*q2.w - q1.y*q2.z + q1.z*q2.y
  result.y =  q1.w*q2.y + q1.x*q2.z + q1.y*q2.w - q1.z*q2.x
  result.z =  q1.w*q2.z - q1.x*q2.y + q1.y*q2.x + q1.z*q2.w
  result.w =  q1.w*q2.w + q1.x*q2.x + q1.y*q2.y + q1.z*q2.z

# Function to rotate vector v using quaternion q
proc transformLeftHand(v: Vector3D, q: Quat): Vector3D =
  let b = vector(q.x, q.y,q.z)
  let b2 = dot(b, b)
  result = scale(v,(q.w*q.w - b2)) + scale(b,(2*dot(v, b))) - scale(cross(b, v),(2*q.w))



#[ Templates ]#






# Make ops that allow mutable compound assignment operation
template makeMutVectorOp(op) =
  proc op*(A: var Vector2D, B: Vector2D) {.inline.} =
    op(A.x, B.x)
    op(A.y, B.y)

  proc op*(A: var Vector3D, B: Vector3D) {.inline.} =
    op(A.x, B.x)
    op(A.y, B.y)
    op(A.z, B.z)

  proc op*(A: var Vector4D, B: Vector4D) {.inline.} =
    op(A.x, B.x)
    op(A.y, B.y)
    op(A.z, B.z)
    op(A.w, B.w)

makeMutVectorOp(`+=`)
makeMutVectorOp(`-=`)
makeMutVectorOp(`*=`)
makeMutVectorOp(`/=`)

# Functions to calculate projection between two vectors (A onto B)
template makeProjs(vector: typedesc) =
  func project*(A, B: vector): vector {.inline.} =
    let scaleVal = dot(A, B)/(getMag(B)^2)
    result = scale(B, scaleVal)

makeProjs(Vector2D)
makeProjs(Vector3D)
makeProjs(Vector4D)

# Functions to calculate the rejection between two vectors (A from B)
template makeRejs(vector: typedesc) =
  func reject*(A, B: vector): vector {.inline.} =
    let magB = getMag(B)
    let dotProd = dot(A, B)
    result = A - scale(B, dotProd/(magB^2))

makeRejs(Vector2D)
makeRejs(Vector3D)
makeRejs(Vector4D)


when isMainModule:
  var v = Vector3D(x: 1.0, y: 0.0, z: 0.0)       # vector along X
  var q = Quat(x: 0.0, y: 0.707, z: 0.0, w: 0.707)  # 90° around Y

  let vRot = transformLeftHand(v, q)
  echo "Original vector: ", v
  echo "Rotated vector:  ", vRot
