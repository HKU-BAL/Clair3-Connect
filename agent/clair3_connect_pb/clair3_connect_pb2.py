"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import runtime_version as _runtime_version
from google.protobuf import symbol_database as _symbol_database
from google.protobuf.internal import builder as _builder
_runtime_version.ValidateProtobufRuntimeVersion(
    _runtime_version.Domain.PUBLIC,
    6,
    31,
    1,
    '',
    'clair3-connect.proto'
)

_sym_db = _symbol_database.Default()

DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n\x14\x63lair3-connect.proto\x12\x11\x63lair3_connect.v1\"\xb1\x01\n\x0ePredictRequest\x12\x12\n\nrequest_id\x18\x01 \x01(\t\x12\x10\n\x08model_id\x18\x02 \x01(\t\x12-\n\x05\x64type\x18\x03 \x01(\x0e\x32\x1e.clair3_connect.v1.TensorDType\x12\r\n\x05shape\x18\x04 \x03(\x03\x12\x0f\n\x07payload\x18\x05 \x01(\x0c\x12\x16\n\x0epayload_sha256\x18\x06 \x01(\t\x12\x12\n\naead_nonce\x18\x07 \x01(\x0c\"\xb3\x02\n\x0fPredictResponse\x12\x12\n\nrequest_id\x18\x01 \x01(\t\x12\x30\n\x06status\x18\x02 \x01(\x0e\x32 .clair3_connect.v1.PredictStatus\x12\x12\n\nerror_code\x18\x03 \x01(\t\x12\x15\n\rerror_message\x18\x04 \x01(\t\x12/\n\x07y_dtype\x18\n \x01(\x0e\x32\x1e.clair3_connect.v1.TensorDType\x12\x0f\n\x07y_shape\x18\x0b \x03(\x03\x12\x11\n\ty_payload\x18\x0c \x01(\x0c\x12\x10\n\x08y_sha256\x18\r \x01(\t\x12\x1c\n\x14input_payload_sha256\x18\x0e \x01(\t\x12\x14\n\x0cmodel_sha256\x18\x0f \x01(\t\x12\x14\n\x0cy_aead_nonce\x18\x10 \x01(\x0c\"\x0f\n\rHealthRequest\"\xb5\x01\n\x0eHealthResponse\x12\x16\n\x0eserver_version\x18\x01 \x01(\t\x12R\n\x11registered_models\x18\x02 \x03(\x0b\x32\x37.clair3_connect.v1.HealthResponse.RegisteredModelsEntry\x1a\x37\n\x15RegisteredModelsEntry\x12\x0b\n\x03key\x18\x01 \x01(\t\x12\r\n\x05value\x18\x02 \x01(\t:\x02\x38\x01*\x8c\x01\n\x0bTensorDType\x12\x1c\n\x18TENSOR_DTYPE_UNSPECIFIED\x10\x00\x12\x15\n\x11TENSOR_DTYPE_INT8\x10\x01\x12\x16\n\x12TENSOR_DTYPE_INT16\x10\x02\x12\x16\n\x12TENSOR_DTYPE_INT32\x10\x03\x12\x18\n\x14TENSOR_DTYPE_FLOAT32\x10\x04*\x98\x01\n\rPredictStatus\x12\x1e\n\x1aPREDICT_STATUS_UNSPECIFIED\x10\x00\x12\x15\n\x11PREDICT_STATUS_OK\x10\x01\x12\x18\n\x14PREDICT_STATUS_RETRY\x10\x02\x12\x19\n\x15PREDICT_STATUS_FAILED\x10\x03\x12\x1b\n\x17PREDICT_STATUS_REJECTED\x10\x04\x32\xb9\x01\n\x16\x43lair3ConnectInference\x12P\n\x07Predict\x12!.clair3_connect.v1.PredictRequest\x1a\".clair3_connect.v1.PredictResponse\x12M\n\x06Health\x12 .clair3_connect.v1.HealthRequest\x1a!.clair3_connect.v1.HealthResponseb\x06proto3')

_globals = globals()
_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, _globals)
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'clair3_connect_pb2', _globals)
if not _descriptor._USE_C_DESCRIPTORS:
  DESCRIPTOR._loaded_options = None
  _globals['_HEALTHRESPONSE_REGISTEREDMODELSENTRY']._loaded_options = None
  _globals['_HEALTHRESPONSE_REGISTEREDMODELSENTRY']._serialized_options = b'8\001'
  _globals['_TENSORDTYPE']._serialized_start=735
  _globals['_TENSORDTYPE']._serialized_end=875
  _globals['_PREDICTSTATUS']._serialized_start=878
  _globals['_PREDICTSTATUS']._serialized_end=1030
  _globals['_PREDICTREQUEST']._serialized_start=44
  _globals['_PREDICTREQUEST']._serialized_end=221
  _globals['_PREDICTRESPONSE']._serialized_start=224
  _globals['_PREDICTRESPONSE']._serialized_end=531
  _globals['_HEALTHREQUEST']._serialized_start=533
  _globals['_HEALTHREQUEST']._serialized_end=548
  _globals['_HEALTHRESPONSE']._serialized_start=551
  _globals['_HEALTHRESPONSE']._serialized_end=732
  _globals['_HEALTHRESPONSE_REGISTEREDMODELSENTRY']._serialized_start=677
  _globals['_HEALTHRESPONSE_REGISTEREDMODELSENTRY']._serialized_end=732
  _globals['_CLAIR3CONNECTINFERENCE']._serialized_start=1033
  _globals['_CLAIR3CONNECTINFERENCE']._serialized_end=1218
