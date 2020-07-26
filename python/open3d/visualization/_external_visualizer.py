import open3d as o3d

__all__ = ['ExternalVisualizer', 'EV']


class ExternalVisualizer:
    """This class allows to send data to an external Visualizer

    Example:
        This example sends a point cloud to the visualizer::

            import open3d as o3d
            import numpy as np
            ev = o3d.visualizer.ExternalVisualizer()
            pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(np.random.rand(100,3)))
            ev.set(pcd)

    Args:
        address: The address where the visualizer is running.
            The default is localhost.
        timeout: The timeout for sending data in milliseconds.
    """

    def __init__(self, address='tcp://127.0.0.1:51454', timeout=10000):
        self.address = address
        self.timeout = timeout

    def set(self,
            obj=None,
            path='',
            time=0,
            layer='',
            objs=None,
            connection=None):
        """Send Open3D objects for visualization to the visualizer.

        Example:
            To quickly send a single object just write::
                ev.set(point_cloud)

            To place the object at a specific location in the scene tree do::
                ev.set(point_cloud, path='group/mypoints', time=42, layer='')
            Note that depending on the visualizer some arguments like time or
            layer may not be supported and will be ignored.

            To pass multiple objects use the ``objs`` keyword argument to pass
            a list::
                ev.set(objs=[point_cloud, mesh, camera])
            Each entry in the list can be a tuple specifying all or some of the
            location parameters::
                ev.set(objs=[(point_cloud,'group/mypoints', 1, 'layer1'),
                             (mesh, 'group/mymesh'),
                             camera
                            ]

        Args:
            obj: A geometry or camera object.

            path: A path describing a location in the scene tree.

            time: An integer time value associated with the object.

            layer: The layer associated with the object.

            objs: List of objects or tuples to pass multiple objects. See the
                example section for usage instructions.

            connection: A connection object to use for sending data. This
                parameter can be used to override the default object.
        """
        if connection is None:
            connection = o3d.utility.Connection(address=self.address,
                                                timeout=self.timeout)
        result = []
        if not obj is None:
            if isinstance(obj, o3d.geometry.PointCloud):
                status = o3d.utility.set_point_cloud(obj,
                                                     path=path,
                                                     time=time,
                                                     layer=layer,
                                                     connection=connection)
                result.append(status)
            elif isinstance(obj, o3d.geometry.TriangleMesh):
                status = o3d.utility.set_triangle_mesh(obj,
                                                       path=path,
                                                       time=time,
                                                       layer=layer,
                                                       connection=connection)
                result.append(status)
            elif isinstance(obj, o3d.camera.PinholeCameraParameters):
                status = o3d.utility.set_legacy_camera(obj,
                                                       path=path,
                                                       time=time,
                                                       layer=layer,
                                                       connection=connection)
                result.append(status)
            else:
                raise Exception("Unsupported object type '{}'".format(
                    str(type(obj))))

        if isinstance(objs, (tuple, list)):
            # item can be just an object or a tuple with path, time, layer, e.g.,
            #   set(objs=[point_cloud, mesh, camera])
            #   set(objs=[(point_cloud,'group/mypoints', 1, 'layer1'),
            #             (mesh, 'group/mymesh'),
            #             camera
            #             ]
            for item in objs:
                if isinstance(item, (tuple, list)):
                    if len(item) in range(1, 5):
                        result.append(self.set(*item, connection=connection))
                else:
                    result.append(self.set(item, connection=connection))

        return all(result)

    def get(self, path='', time=0, layer='', connection=None):
        """Get an object from the visualizer.


        Args:
            path: A path describing a location in the scene tree.

            time: An integer time value associated with the object.

            layer: The layer associated with the object.

            connection: A connection object to use for sending data. This
                parameter can be used to override the default object.
        """
        if connection is None:
            connection = o3d.utility.Connection(address=self.address,
                                                timeout=self.timeout)
        obj = o3d.utility.get_mesh_data(path, time, layer)
        return obj

    def set_time(self, time):
        """Sets the time in the external visualizer

        Args:
            time: The time value
        """
        connection = o3d.utility.Connection(address=self.address,
                                            timeout=self.timeout)
        return o3d.utility.set_time(time, connection)

    def set_active_camera(self, path):
        """Sets the time in the external visualizer

        Args:
            path: A path describing a location in the scene tree.
        """
        connection = o3d.utility.Connection(address=self.address,
                                            timeout=self.timeout)
        return o3d.utility.set_active_camera(path, connection)


# convenience default external visualizer
EV = ExternalVisualizer()
