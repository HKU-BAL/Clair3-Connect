use anyhow::Result;
use tch::Tensor;

use crate::torchscript::LoadedModel;

pub struct Worker {

    pub model: LoadedModel,
}

impl Worker {

    pub fn new(model: LoadedModel) -> Self {
        Self { model }
    }

    pub fn forward_batch(&self, x: &Tensor) -> Result<Tensor> {
        self.model.forward(x)
    }
}
