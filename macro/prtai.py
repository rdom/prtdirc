import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import numpy as np

import tensorflow as tf
print("TensorFlow version:", tf.__version__)
from tensorflow.keras.layers import Dense, Flatten, Conv2D, Discretization
from tensorflow.keras import Model
from array import array

import ROOT
from ROOT import TFile, TTree, gROOT, addressof

stat = int(sys.argv[1])
        
def load_data(path="data.npz"):
    
    with np.load(path, allow_pickle=True) as f:
        x_train, y_train = f["x_train"], f["y_train"]
        x_test, y_test = f["x_test"], f["y_test"]
        
        return (x_train, y_train), (x_test, y_test)


(x_train, y_train), (x_test, y_test) = load_data("data_stat_16x32_t_40000.npz")
# x_train, x_test = x_train / 30.0, x_test / 30.0

x_train = x_train[:stat]
y_train = y_train[:stat]

# assert x_train.shape == (20000, 8,64)
# assert x_test.shape == (2000, 8,64)
# assert y_train.shape == (20000,1)
# assert y_test.shape == (2000,1)
    
# Add a channels dimension
x_train = x_train[..., tf.newaxis].astype("float")
x_test = x_test[..., tf.newaxis].astype("float")


# tf.print(x_train, summarize=-1)

train_ds = tf.data.Dataset.from_tensor_slices((x_train, y_train)).shuffle(100).batch(32)
test_ds = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(32)

class PrtNN(Model):
  def __init__(self):
    super(PrtNN, self).__init__()
    # self.conv1 = Conv2D(2, 4, activation='relu')
    self.disc = Discretization(bin_boundaries=[0.001])
    self.flatten = Flatten()
    # self.d1 = Dense(5, input_shape = 512, activation='relu')
    self.d2 = Dense(5)

  def call(self, x):
      # x = self.conv1(x)
      x = self.disc(x)
      x = self.flatten(x)
      # x = self.d1(x)
      return self.d2(x)

  def model(self):
      x = tf.keras.Input(shape=(16,32,1))
      return Model(inputs=[x], outputs=self.call(x))
  

model = PrtNN()
loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
train_loss = tf.keras.metrics.Mean(name='train_loss')
train_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='train_accuracy')
test_loss = tf.keras.metrics.Mean(name='test_loss')
test_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='test_accuracy')

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

tf.keras.utils.plot_model(model, to_file='model.png', show_shapes=True)


model.model().summary()
model.save('models/prtai')


print("Accuracy = ", faccuracy)

rf = TFile("res_" + str(stat)  + ".root", "RECREATE")
tree = TTree("T", "prtai")

eff = array( 'f', [ 0 ] )
sta = array( 'i', [ 0 ] )



tree.Branch( 'sta', sta, 'sta/I' )
tree.Branch( 'eff', eff, 'eff/F' )

eff[0] = float(faccuracy)
sta[0] = stat

tree.Fill()
tree.Print()
tree.Write()



