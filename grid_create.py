import gmsh
import sys
import os
import subprocess

def generate_mesh():
    # Инициализация gmsh
    gmsh.initialize()
    gmsh.model.add("model")

    # Параметры сетки
    d1 = 50
    d2 = 60
    d3 = 100
    d4 = 60

    LX = 0.02
    dy = 0.002
    LY = 0.01 + dy

    p1 = LX / d1
    p2 = LY / d2
    p3 = LY / d3
    p4 = LY / d4

    # Создаем точки
    point1 = gmsh.model.geo.addPoint(0, 0, 0, p1)
    point2 = gmsh.model.geo.addPoint(0, dy, 0, p2)
    point3 = gmsh.model.geo.addPoint(0, LY, 0, p3)
    point4 = gmsh.model.geo.addPoint(LX, LY, 0, p3)
    point5 = gmsh.model.geo.addPoint(LX, dy, 0, p2)
    point6 = gmsh.model.geo.addPoint(LX, 0, 0, p1)
    point21 = gmsh.model.geo.addPoint(0, LY - LY / 5, 0, p4)
    point41 = gmsh.model.geo.addPoint(LX, LY - LY / 5, 0, p4)
    point32 = gmsh.model.geo.addPoint(LX - LX / 4, LY, 0, p4)
    point31 = gmsh.model.geo.addPoint(LX / 4, LY, 0, p4)

    # Создаем линии и замкнутые линейные петли
    line1 = gmsh.model.geo.addLine(point1, point2)
    line2 = gmsh.model.geo.addLine(point2, point21)
    line3 = gmsh.model.geo.addLine(point21, point3)
    line4 = gmsh.model.geo.addLine(point3, point31)
    line5 = gmsh.model.geo.addLine(point31, point32)
    line6 = gmsh.model.geo.addLine(point32, point4)
    line7 = gmsh.model.geo.addLine(point4, point41)
    line8 = gmsh.model.geo.addLine(point41, point5)
    line9 = gmsh.model.geo.addLine(point5, point6)
    line10 = gmsh.model.geo.addLine(point6, point1)
    line11 = gmsh.model.geo.addLine(point5, point2)

    loop1 = gmsh.model.geo.addCurveLoop([line2, line3, line4, line5, line6, line7, line8, line11])
    loop2 = gmsh.model.geo.addCurveLoop([line1, -line11, line9, line10])
    surf1 = gmsh.model.geo.addPlaneSurface([loop1])
    surf2 = gmsh.model.geo.addPlaneSurface([loop2])

    # Синхронизация и создание физических групп
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [surf1], 1)
    gmsh.model.addPhysicalGroup(2, [surf2], 2)
    gmsh.model.addPhysicalGroup(1, [line4, line5, line6], 3)
    gmsh.model.addPhysicalGroup(1, [line1, line2, line3], 4)
    gmsh.model.addPhysicalGroup(1, [line7, line8, line9], 5)
    gmsh.model.addPhysicalGroup(1, [line10], 6)

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    output_dir_model = "./P_model/"
    output_dir_fen = "./P_Fen/"
    os.makedirs(output_dir_model, exist_ok=True)
    os.makedirs(output_dir_fen, exist_ok=True)

    mesh_file_base = "mesh"
    gmsh.write(f"{output_dir_model}/{mesh_file_base}.msh")

    try:
        subprocess.run(['dolfin-convert', f'{output_dir_model}/{mesh_file_base}.msh', f'{output_dir_fen}/{mesh_file_base}.xml'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Ошибка конвертации: {e}")

    gmsh.finalize()

if __name__ == "__main__":
    generate_mesh()
