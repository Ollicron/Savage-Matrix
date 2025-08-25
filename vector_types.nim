#[
The purpose of this module is to define types used by the graphics vectors module and define operators to access elements them. It's essential that the vectors are objects placed on the stack as opposed to the heap because performance is more essential than feasability. It's easier to have all the functions laid out here to reduce compilation time.
]#

type
  Vector2D* = object
    x*, y*: float

  Vector3D* = object of RootObj
    x*, y*, z*: float

  Vector4D* = object
    x*, y*, z*, w*: float

  Point3D* = object of Vector3D


func vector*(A: int|float, B: int|float): Vector2D =
  result = Vector2D(x: float(A), y: float(B))

func vector*(A: int|float, B: int|float, C: int|float): Vector3D =
  result = Vector3D(x: float(A), y: float(B), z: float(C))

func vector*(A: int|float, B: int|float, C: int|float, D: int|float): Vector4D =
  result = Vector4D(x: float(A), y: float(B), z: float(C), w: float(D))

# New array-based constructors
func vector*(arr: array[2, SomeNumber]): Vector2D =
  Vector2D(x: float(arr[0]), y: float(arr[1]))

func vector*(arr: array[3, SomeNumber]): Vector3D =
  Vector3D(x: float(arr[0]), y: float(arr[1]), z: float(arr[2]))

func vector*(arr: array[4, SomeNumber]): Vector4D =
  Vector4D(x: float(arr[0]), y: float(arr[1]), z: float(arr[2]), w: float(arr[3]))

# Constructors for Point3D
func point3D*(A: int|float, B: int|float, C: int|float): Point3D =
  result = Point3D(x: float(A), y: float(B), z: float(C))

func point3D*(vector:Vector3D): Point3D =
  result = Point3D(x: vector.x, y: vector.y, z: vector.z)


func `[]`*(vector: Vector2D, i: int): float {.discardable inline.} =
  case i
  of 0: result = vector.x
  of 1: result = vector.y
  else: discard

func `[]`*(vector: Vector3D, i: int): float {.discardable inline.} =
  case i
  of 0: result = vector.x
  of 1: result = vector.y
  of 2: result = vector.z
  else: discard

func `[]`*(vector: Vector4D, i: int): float {.discardable inline.} =
  case i
  of 0: result = vector.x
  of 1: result = vector.y
  of 2: result = vector.z
  of 3: result = vector.w
  else: discard

proc `[]=`*(vector: var Vector2D, i: float|int) {.inline.} =
  case i
  of 0: vector.x = float(i)
  of 1: vector.y = float(i)
  else: discard

proc `[]=`*(vector: var Vector3D, i: float|int) {.inline.} =
  case i
  of 0: vector.x = float(i)
  of 1: vector.y = float(i)
  of 2: vector.z = float(i)
  else: discard

proc `[]=`*(vector: var Vector4D, i: float|int){.inline.} =
  case i
  of 0: vector.x = float(i)
  of 1: vector.y = float(i)
  of 2: vector.z = float(i)
  of 3: vector.w = float(i)
  else: discard
