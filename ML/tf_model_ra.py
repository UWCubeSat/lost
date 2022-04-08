import tensorflow as tf
import pandas as pd
import tensorflow.keras as keras
from constants import *
import wandb
from wandb.keras import WandbCallback
wandb.init(project="ra-proj", entity="sathvikc")
wandb.config = {
  "learning_rate": LEARNING_RATE,
  "epochs": NUM_EPOCHS,
  "batch_size": BATCH_SIZE,
  "num_images": NUM_IMAGES
}
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

data = pd.read_csv(PATH_TO_DATA_CSV)

def decode_img(img):
  # Convert the compressed string to a 3D uint8 tensor
  img = tf.io.decode_png(img, channels=3)
  # Resize the image to the desired size
  return tf.image.resize(img, [1000, 1000])


all_image_paths = data['image_path'].to_numpy()
all_image_paths = DATA_PATH + '/' + all_image_paths
all_image_labels = data['ra'].to_numpy()

def load_and_preprocess_image(path):
  image_string = tf.compat.as_str_any(path)
  image_string = tf.io.read_file(path)
  img = tf.io.decode_png(image_string, channels=3)
  return tf.image.resize(img, [1000, 1000])

def load_and_preprocess_from_path_labels(path, label):
  return load_and_preprocess_image(path), label

ds = tf.data.Dataset.from_tensor_slices((all_image_paths, all_image_labels))
star_csv_dataset = ds.map(load_and_preprocess_from_path_labels, num_parallel_calls=tf.data.AUTOTUNE)

def configure_for_performance(ds):
  ds = ds.cache()
  ds = ds.shuffle(buffer_size=1000)
  ds = ds.batch(BATCH_SIZE)
  ds = ds.prefetch(buffer_size=tf.data.AUTOTUNE)
  return ds

# star_csv_dataset = configure_for_performance(star_csv_dataset)
train_size = int(0.7 * NUM_IMAGES)
val_size = int(0.15 * NUM_IMAGES)
test_size = int(0.15 * NUM_IMAGES)
star_csv_dataset.shuffle(buffer_size = 1000, reshuffle_each_iteration=False)
train_dataset = star_csv_dataset.take(train_size)
test_dataset = star_csv_dataset.skip(train_size)
val_dataset = test_dataset.skip(val_size)
test_dataset = test_dataset.take(test_size)

train_dataset = configure_for_performance(train_dataset)
val_dataset = configure_for_performance(val_dataset)
test_dataset = configure_for_performance(test_dataset)

callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5)
model = keras.Sequential([
  tf.keras.layers.Conv2D(8, 3, activation='relu'),
  tf.keras.layers.MaxPooling2D(),
  tf.keras.layers.Conv2D(8, 3, activation='relu'),
  tf.keras.layers.MaxPooling2D(),
  tf.keras.layers.Flatten(),
  tf.keras.layers.Dense(1, activation = 'sigmoid'), 
  tf.keras.layers.Lambda(lambda x : x * 365)
])

model.compile(optimizer='adam',
    loss = tf.keras.losses.MeanAbsoluteError(), 
    metrics = [
      tf.keras.metrics.MeanAbsoluteError(), 
      tf.keras.metrics.MeanSquaredError(), 
      tf.keras.metrics.MeanAbsolutePercentageError()])

for image, label in train_dataset.take(1):
  print(model.predict(image))
  print(label)
  mse = tf.keras.losses.MeanAbsoluteError()
  print(mse(label, model.predict(image)))

model.fit(train_dataset, epochs = NUM_EPOCHS, validation_data = val_dataset, callbacks=[WandbCallback(), callback])
model.evaluate(test_dataset, callbacks=[WandbCallback()])
