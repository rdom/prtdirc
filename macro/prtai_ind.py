import sys
import os
import re
import numpy as np
from array import array
from ROOT import TFile, TTree

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
from tensorflow.keras import Model
from tensorflow.keras.layers import Dense, Flatten, Conv2D, Discretization, Normalization, Dropout

# print("TensorFlow version:", tf.__version__)

infile = sys.argv[1]
   
def load_data(path="data.npz"):    
    with np.load(path, allow_pickle=True) as f:
        data_x, data_y = f["x_train"].astype("int32"), f["y_train"]
        
        return (data_x, data_y)


# (data_x, data_y) = load_data("ind_97_30000.npz")
(data_x, data_y) = load_data(infile)


if(len(sys.argv) > 2):
    stat = int(sys.argv[2])
else :
    stat = data_x.shape[0]
       
print("Number of events: ", stat)
data_x = data_x[:stat]
data_y = data_y[:stat]
                                                    

# # Add a channels dimension
# data_x = data_x[..., tf.newaxis].astype("int32")
# x_test = x_test[..., tf.newaxis].astype("int32")

# with np.printoptions(precision=0, linewidth=300, edgeitems=100):
#     print(data_x)

# exit()

nbatches = 64
# data_ds = tf.data.Dataset.from_tensor_slices((data_x, data_y)).shuffle(stat).batch(32)
data_ds = tf.data.Dataset.from_tensor_slices((data_x, data_y)).batch(nbatches)
train_ds, test_ds = tf.keras.utils.split_dataset(data_ds, left_size=0.75)

print("Training on:", nbatches * int(train_ds.cardinality()),
      " Testing on: ", nbatches * int(test_ds.cardinality()))


# for x, y in train_ds:
#     print(x, y)

# tf.config.run_functions_eagerly(True)

timebins = 50

class PrtNN(Model):
    def __init__(self):
        super(PrtNN, self).__init__()
        self.conv1 = Conv2D(filters=4,kernel_size=(2,2),
                            kernel_initializer='glorot_uniform',
                            activation='relu')
        self.maxpool = tf.keras.layers.MaxPooling2D(pool_size=(16, 5))
        self.disc = Discretization(bin_boundaries=[0.001])
        self.norm = Normalization(axis=None)
        self.flatten = Flatten()
        self.d1 = Dense(5, activation='relu')
        self.d2 = Dense(5)
        self.drop = Dropout(0.2)

    def call(self, x):

        batches = x.shape[0]
        if batches is None:
            batches = 1

        # with np.printoptions(precision=0, linewidth=300, edgeitems=100):
        #     print(x)
        # exit()
                
        # convert indexes into tensor
        b = tf.range(batches)               
        b = tf.expand_dims(b, -1)
        b = tf.repeat(b, repeats=100, axis=-1)
        b = tf.expand_dims(b, -1)
        x = tf.concat((b, x), axis=-1)          
        # with np.printoptions(precision=0, linewidth=300, edgeitems=100):
        #     print(x)

        # new way
        whits = tf.ones([batches,98], tf.float32)
        wtracks = tf.fill([batches,2], 2.0)        
        xhits, xtracks = tf.split(x, [98, 2], 1)        
        zhits = tf.zeros([batches,512,timebins], tf.float32)
        ztracks = tf.zeros([batches,20,timebins], tf.float32)
        xhits = tf.tensor_scatter_nd_update(zhits, xhits, whits)
        xtracks = tf.tensor_scatter_nd_update(ztracks, xtracks, wtracks)
        x = tf.concat((xhits, xtracks), axis=1) 
        # with np.printoptions(precision=0, linewidth=300, edgeitems=100):
        #     print(x)

        # old way
        # ones = tf.ones([batches,100], tf.float32)
        # z = tf.zeros([batches,512 + 20,timebins], tf.float32) # 512,50
        # x = tf.tensor_scatter_nd_update(zhit, x, ones)


        # x = tf.expand_dims(x, -1)
        
        # x = self.norm(x)
        # x = self.conv1(x)
        # x = self.d1(x)
        # x = self.maxpool(x)
        # x = self.disc(x)
        x = self.flatten(x)
        # x = self.drop(x)
        # x = self.d1(x)
        return self.d2(x)

    def model(self):
        # x = tf.keras.Input(shape=(16,32,100))
        x = tf.keras.Input(shape=(100,3))
        x = tf.cast(x, tf.int32)
        return Model(inputs=[x], outputs=self.call(x))
    

model = PrtNN()
loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
train_loss = tf.keras.metrics.Mean(name='train_loss')
train_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='train_accuracy')
test_loss = tf.keras.metrics.Mean(name='test_loss')
test_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='test_accuracy')

model.compile(optimizer='adam',
              loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
              metrics=[tf.keras.metrics.SparseCategoricalAccuracy()])


@tf.function
def train_step(images, labels):
    with tf.GradientTape() as tape:
        predictions = model(images, training=True)
        loss = loss_object(labels, predictions)
        
    gradients = tape.gradient(loss, model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))
    
    train_loss(loss)
    train_accuracy(labels, predictions)
    
@tf.function
def test_step(images, labels):
    predictions = model(images, training=False)
    t_loss = loss_object(labels, predictions)
    
    test_loss(t_loss)
    test_accuracy(labels, predictions)

  
faccuracy = 0
EPOCHS = 10

for epoch in range(EPOCHS):
    # Reset the metrics at the start of the next epoch
    train_loss.reset_states()
    train_accuracy.reset_states()
    test_loss.reset_states()
    test_accuracy.reset_states()
    
    for images, labels in train_ds:
        train_step(images, labels)

    for test_images, test_labels in test_ds:
        test_step(test_images, test_labels)
            
    faccuracy = test_accuracy.result() * 100
    print(
        f'Epoch {epoch + 1}, '
        f'Loss: {train_loss.result()}, '
        f'Accuracy: {train_accuracy.result() * 100}, '
        f'Test Loss: {test_loss.result()}, '
        f'Test Accuracy: {test_accuracy.result() * 100}'
    )

# tf.keras.utils.plot_model(model, to_file='model.png', show_shapes=True)

model.model().summary()
model.save('models/' + infile)


print("Accuracy = ", faccuracy)

outfile = infile.replace("npz","ires.root")
angle = float(re.findall(r'\d+.\d+', infile)[0])
print(outfile, " theta = ", angle, " eff = ", float(faccuracy))

# outfile = "res_" + str(stat)  + ".root"
rf = TFile(outfile, "RECREATE")
tree = TTree("T", "prtai")

eff = array( 'f', [ 0 ] )
sta = array( 'i', [ 0 ] )
theta = array( 'f', [ 0 ] )

tree.Branch( 'sta', sta, 'sta/I' )
tree.Branch( 'eff', eff, 'eff/F' )
tree.Branch( 'theta', theta, 'theta/F' )

eff[0] = float(faccuracy)
sta[0] = stat
theta[0] = angle

tree.Fill()
tree.Write()



