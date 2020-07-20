# ParFlow Diagnostics

This Diagnositics package provides functions to calculated a global or local mass balance based on ParFlow output

**Diagnostics.py:**

SubsurfaceStorage(self, Press, Satur) -> subsurface_storage

VolumetricMoisture(self, Satur) -> volumetric_moisture´

TopLayerPressure(self, Press) -> top_layer_pressure

SurfaceStorage(self,Toplayerpress) -> surface_storage

OverlandFlow(self,Toplayerpress) -> oflowx, oflowy

NetLateralOverlandFlow(self, overland_flow_x, overland_flow_y) -> net_overland_flow

SubsurfaceFlow(self, Press, Krel) -> flowleft, flowright, flowfront, flowback, flowbottom, flowtop

VanGenuchten(self,Press) -> satur, krel