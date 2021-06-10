from typing import List, Dict

import numpy as np

from src.models.annotation_model import AnnotationModel
from src.util.feature import Feature


class TransitionModel:
    def __init__(self, annotation_models: List[AnnotationModel]):
        self._models = annotation_models
        self.feature_slls: Dict[str, List[Feature]] = AnnotationModel.merge(annotation_models)
        self.transition_model: np.ndarray = self._calculate_transition_model()

    def _calculate_transition_model(self) -> np.ndarray:
        pass
