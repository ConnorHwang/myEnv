import bpy
import glob
import math
import numpy as np

pi = math.pi
   
filename=glob.glob("/p/lustre1/hwang8/dfdstl_case1/*.stl")
filename.sort()
min=0
max=len(filename)-1

counter=min # Counter to count the number of STL files
#frame = 1 # the actual frame number in the blender (your animation)
frame_max = 3 # max(frame) <= max(counter)

savepath = "/p/lustre1/hwang8/rendering_case0/"
savename = "case7"
#------ LOOP   SETTING -------
index_init = 0
index_final = 341
frame = index_init + 1
#------ CAMERA SETTING -------
initCamLocX=10
initCamLocY=0
initCamLocZ=-0.6
initCamRotX=-90
initCamRotY=0
initCamRotZ=-90
myLightX=6
myLightY=5
myLightZ=-7
camFocalLength=25
# Test for loop
cameraBegin=initCamLocX
cameraEnd=cameraBegin+5
#Nc=index_final-index_init
Nc=341
cameraPos=np.linspace(cameraBegin,cameraEnd,Nc)
#-----------------------------

#notdone = True

bpy.ops.object.delete()

# Remove the default light
bpy.ops.object.select_by_type(type='LIGHT')
bpy.ops.object.delete(use_global=False)
light_data = bpy.data.lights.new(name="my_light_data",type='POINT')
light_object = bpy.data.objects.new(name="my_light",object_data=light_data)
bpy.context.collection.objects.link(light_object)
light_object.location=(myLightX,myLightY,myLightZ)
light_data=light_object.data
light_data.energy=2500
light_data.shadow_soft_size=10

# Set the camera
# Remove the default camera
bpy.ops.object.select_by_type(type='CAMERA')
bpy.ops.object.delete(use_global=False)
camera_data = bpy.data.cameras.new(name="my_camera_data")
camera_object = bpy.data.objects.new("my_camera",camera_data)
camera_data.lens=camFocalLength
bpy.context.collection.objects.link(camera_object)
camera_object.location=(initCamLocX,initCamLocY,initCamLocZ)
camera_object.rotation_mode = 'XYZ'
camera_object.rotation_euler[0] = initCamRotX * pi / 180
camera_object.rotation_euler[1] = initCamRotY * pi / 180
camera_object.rotation_euler[2] = initCamRotZ * pi / 180

# Add planes
## 1st plane
bpy.ops.mesh.primitive_plane_add(size=30,enter_editmode=False,location=(0,0,-20))
bpy.context.active_object.name='plane_bottom'
bpy.ops.material.new()
bpy.data.materials[1].name='plane_color'
bpy.data.materials['plane_color'].diffuse_color=(0.2,0.2,0.2,1)
bpy.data.materials['plane_color'].node_tree.nodes["Principled BSDF"].inputs[0].default_value = (0.1, 0.1, 0.1, 1.0)


bpy.context.object.data.materials.append(bpy.data.materials['plane_color'])
# 2nd plane
bpy.ops.mesh.primitive_plane_add(size=100,enter_editmode=False,location=(0,-15,0))
plane = bpy.context.active_object
plane.name='plane_right'
plane.rotation_mode = 'XYZ'
plane.rotation_euler[0] = 90*pi/180
plane.rotation_euler[1] = 0
plane.rotation_euler[2] = 0
bpy.context.object.data.materials.append(bpy.data.materials['plane_color'])
# 3rd plane
bpy.ops.mesh.primitive_plane_add(size=100,enter_editmode=False,location=(0,15,0))
plane = bpy.context.active_object
plane.name='plane_left'
plane.rotation_mode = 'XYZ'
plane.rotation_euler[0] = -90*pi/180
plane.rotation_euler[1] = 0
plane.rotation_euler[2] = 0
bpy.context.object.data.materials.append(bpy.data.materials['plane_color'])
# 4th plane
bpy.ops.mesh.primitive_plane_add(size=100,enter_editmode=False,location=(-1,0,0))
plane = bpy.context.active_object
plane.name='plane_nozzle'
plane.rotation_mode = 'XYZ'
plane.rotation_euler[0] = 0
plane.rotation_euler[1] = 90*pi/180
plane.rotation_euler[2] = 0
bpy.context.object.data.materials.append(bpy.data.materials['plane_color'])
# 5th plane
bpy.ops.mesh.primitive_plane_add(size=100,enter_editmode=False,location=(50,0,0))
plane = bpy.context.active_object
plane.name='plane_far'
plane.rotation_mode = 'XYZ'
plane.rotation_euler[0] = 0
plane.rotation_euler[1] = 90*pi/180
plane.rotation_euler[2] = 0
bpy.context.object.data.materials.append(bpy.data.materials['plane_color'])
#-----------

# Add Cylinder
bpy.ops.mesh.primitive_cylinder_add(vertices=128,radius=0.4,depth=1,location=(-0.5,0,0),rotation=(0,90*pi/180,0),scale=(1,1,1))
cylinder = bpy.context.active_object
cylinder.name='nozzle'
bpy.ops.material.new()
bpy.data.materials[2].name = 'CylColor'
bpy.data.materials['CylColor'].diffuse_color=(0.2,0.2,0.2,1)
bpy.data.materials['CylColor'].node_tree.nodes["Principled BSDF"].inputs[0].default_value = (0.0, 0.0, 0.0, 0.5)
bpy.data.materials['CylColor'].node_tree.nodes["Principled BSDF"].inputs[19].default_value = 1.0
#bpy.data.materials['CylColor'].node_tree.nodes["Principled BSDF"].inputs[7].default_value = 1.0
bpy.context.object.data.materials.append(bpy.data.materials['CylColor'])
#-----------

# Change rendering environment t0 CYCLE
bpy.context.scene.render.engine='CYCLES'
#bpy.context.scene.cycles.samples=256
 
#----------------
# Test sigle files
#---------------- 
#bpy.ops.object.select_pattern(pattern="n6jet*")
#bpy.ops.import_mesh.stl(filepath=filename[counter])
#jet_surface = bpy.context.active_object
##jet_surface.location=(0.0,-1.65,-1.65)
#bpy.ops.transform.resize(value=(100,100,100))
#bpy.ops.object.shade_smooth()

##bpy.data.materials["Material"]
#material = bpy.data.materials.new(name="liquid_fuel")
#material.use_nodes = True
#material.node_tree.nodes.new(type='ShaderNodeBsdfGlass')
#inp = material.node_tree.nodes['Material Output'].inputs['Surface']
#outp = material.node_tree.nodes['Glass BSDF'].outputs['BSDF']
#material.node_tree.links.new(inp,outp)
#jet_surface.data.materials.append(material)

#----------------
# Loop
#----------------
for ii in range(index_init,index_final):
    
    # Import STL file
    bpy.ops.object.select_pattern(pattern="case*")
    bpy.ops.import_mesh.stl(filepath=filename[ii]) # ii instead of counter
    jet_surface = bpy.context.active_object
#    jet_surface.location=(0.0,-1.65,-1.65)
    bpy.ops.transform.resize(value=(100,100,100))
    bpy.ops.object.shade_smooth()

#    material = bpy.data.materials.new(name="liquid_fuel")
#    material.use_nodes = True
#    material.node_tree.nodes.new(type='ShaderNodeBsdfGlass')
#    inp = material.node_tree.nodes['Material Output'].inputs['Surface']
#    outp = material.node_tree.nodes['Glass BSDF'].outputs['BSDF']
#    material.node_tree.links.new(inp,outp)
#    jet_surface.data.materials.append(material)
    
    # Principle BSDF
    bpy.ops.material.new()
    bpy.data.materials[3].name = 'fuel'
    bpy.data.materials['fuel'].diffuse_color=(0.2,0.2,0.2,1)
    bpy.data.materials['fuel'].node_tree.nodes["Principled BSDF"].inputs[0].default_value = (0.0, 0.2, 1.0, 1.0)
    bpy.data.materials['fuel'].node_tree.nodes["Principled BSDF"].inputs[7].default_value = 0.5 # roughness
    bpy.data.materials['fuel'].node_tree.nodes["Principled BSDF"].inputs[15].default_value = 0.0 # transmission
    bpy.data.materials['fuel'].node_tree.nodes["Principled BSDF"].inputs[19].default_value = 1.0 # alpha
    bpy.context.object.data.materials.append(bpy.data.materials['fuel'])

    # File path
    if frame < 10 :
        filepath=savepath+savename+"."+"000"+str(frame)+".png"
    elif frame < 100 :
        filepath=savepath+savename+"."+"00"+str(frame)+".png"  
    elif frame < 1000 :
        filepath=savepath+savename+"."+"0"+str(frame)+".png"  
    else :
        filepath=savepath+savename+"."+str(frame)+".png"  
        
    # New camera location
    camera_object.location[0]=cameraPos[ii]
    
    # Render
    bpy.context.scene.camera=camera_object
    scene = bpy.context.scene
    bpy.ops.render.render()
    bpy.data.images['Render Result'].save_render(filepath)

    # Increase counts
    frame = frame + 1
    counter = counter + 1
    
    print("finish step=",frame,counter)
    print('Current i: '+str(ii)+', range: '+str(index_init)+'~'+str(index_final))
    
    # Delete jet surface
    bpy.ops.object.delete()

#----------------    
# End of the loop  
