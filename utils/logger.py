import tensorflow as tf


# class Logger(object):

#     def __init__(self, log_dir):
#         """Create a summary writer logging to log_dir."""
#         #self.writer = tf.compat.v1.summary.FileWriter(log_dir)
#         self.writer = tf.summary.create_file_writer(log_dir)


#     def scalar_summary(self, tag, value, step):
#         """Log a scalar variable."""
#         summary = tf.compat.v1.Summary(value=[tf.Summary.Value(tag=tag, simple_value=value)])
#         self.writer.add_summary(summary, step)



class Logger:
    def __init__(self, log_dir):
        self.writer = tf.summary.create_file_writer(log_dir)

    def scalar_summary(self, tag, value, step):
        with self.writer.as_default():
            tf.summary.scalar(name=tag, data=value, step=step)
            self.writer.flush()