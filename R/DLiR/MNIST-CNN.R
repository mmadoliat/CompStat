library(keras3)

# Model / data parameters
num_classes <- 10
input_shape <- c(28, 28, 1)

# Load the data and split it between train and test sets
c(c(x_train, y_train), c(x_test, y_test)) %<-% dataset_mnist()

# Scale images to the [0, 1] range
x_train <- x_train / 255
x_test <- x_test / 255

# Make sure images have shape (28, 28, 1)
x_train <- op_expand_dims(x_train, -1)
x_test <- op_expand_dims(x_test, -1)

# convert class vectors to binary class matrices
y_train <- to_categorical(y_train, num_classes)
y_test <- to_categorical(y_test, num_classes)

# Define Model -----------------------------------------------------------

# Define model
model <- keras_model_sequential(input_shape = input_shape)
model |>
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu") |>
  layer_max_pooling_2d(pool_size = c(2, 2)) |>
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") |>
  layer_max_pooling_2d(pool_size = c(2, 2)) |>
  layer_flatten() |>
  layer_dropout(rate = 0.5) |>
  layer_dense(units = num_classes, activation = "softmax")

# Compile model
batch_size <- 128
epochs <- 15

model |> compile(
  loss = "categorical_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

# Train model
model |> fit(
  x_train, y_train,
  batch_size = batch_size,
  epochs = epochs,
  validation_split = 0.1
)




scores <- model %>% evaluate(
  x_test, y_test, verbose = 0
)

# Output metrics
score <- model |> evaluate(x_test, y_test, verbose = 0)
score

