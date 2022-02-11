resource "aws_sfn_state_machine" "sfn_state_machine" {
  name       = var.sfn_state_machine_name 
  role_arn   = aws_iam_role.iam_for_sfn.arn

  definition = templatefile(var.sfn_statemachine_json_file)
}
